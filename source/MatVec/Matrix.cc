#include "Matrix.hh"
#include "Vector.hh"
#include "opdefs.hh"

#include <fstream>
#include <string>
#include <def_build_type_options.hh>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "Utils/boost-serialization.hh"
#include "Utils/tools.hh"

#include "BLASLAPACKInterface.hh"

using boost::tokenizer;

namespace CoupledField
{      

  template<class TYPE>
  Matrix<TYPE>::Matrix () :
    DenseMatrix(),
    size_row_(0),
    size_col_(0),
    data_(NULL)
  { }


  template<class TYPE>
  Matrix<TYPE>::Matrix (const UInt nRows, const UInt nCols) :
    DenseMatrix(),
    size_row_(nRows),
    size_col_(nCols),
    data_(new TYPE* [size_row_])
  {
#ifdef CHECK_INDEX 
    if (nRows <= 0 || nCols <= 0)
      EXCEPTION("invalid dimension");
#endif

    data_[0]=new TYPE[size_col_*size_row_];

    for (UInt k=1; k < size_row_; k++) 
      data_[k]=data_[k-1]+size_col_;

    Init();
  }


  template<class TYPE>
  Matrix<TYPE>::Matrix (const UInt nRows,const Vector<TYPE> * const x) :
    DenseMatrix(),
    size_row_(nRows),
    size_col_(x[0].size_),
    data_(new TYPE* [size_row_])
  {
#ifdef CHECK_INDEX
    if (nRows <= 0) EXCEPTION("invalid dimension");
    if (size_col_ == 0) EXCEPTION("invalid dimension");
#endif 

    UInt k,kk;
#ifdef CHECK_INDEX
    for (k=1; k < size_row_; k++)
    
      { if (x[k].size_!=size_col_)  
          EXCEPTION(" Not all vectors for initialization have the same size" );
      }
#endif
  
    for (k=0; k < size_row_; k++)
      for (kk=0; kk<size_col_; kk++)
        data_[k][kk]=x[k][kk];
  }

  template<class TYPE>
  Matrix<TYPE>::Matrix (const Matrix<TYPE> &x) :
    DenseMatrix(),
    size_row_(x.size_row_),
    size_col_(x.size_col_),
    data_(NULL)
  {
    // We are able to copy empty matrices!
    if (size_row_ > 0 && size_col_ > 0 )
    {
      data_ = new TYPE * [size_row_];
      data_[0]=new TYPE[size_row_ * size_col_];
    }
 
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k]=x.data_[0][k];
    
    for(UInt k = 1; k < size_row_; ++k)
      data_[k]=data_[k-1]+size_col_;        
  }

  template<class TYPE>
  Matrix<TYPE>::~Matrix ()
  {
    if (data_ != NULL)
    {
      delete[] data_[0];
      delete[] data_;
      data_= NULL;
    }
  }


  template<class TYPE>
  std::string Matrix<TYPE>::ToXMLFormat(const std::string& name, int n_offset) const
  {
    std::string offset(std::max(n_offset,0), ' ');

    std::ostringstream os;

    bool is_complex = boost::is_same<TYPE, std::complex<double> >::value;

    os << "<" << name << " dim1=\"" << size_row_ << "\" dim2=\"" << size_col_ << "\">";
    os << std::endl << offset << "  " << (is_complex ? "<complex>" : "<real>");

    for(unsigned int r = 0; r < size_row_; ++r)
    {
      os << std::endl << offset << "    ";
      for(unsigned int c = 0; c < size_col_; ++c)
      {
        os << std::scientific;
        os.precision(6);
        os.width(13);
        os << (IsNoise(data_[r][c]) ? 0.0 : data_[r][c]);
        if(c < size_col_ - 1) os << " ";
      }
    }
    os << std::endl << offset << "  " << (is_complex ? "</complex>" : "</real>");
    os << std::endl << offset << "</" << name << ">";

    return os.str();
  }

  template<class TYPE>
  std::string Matrix<TYPE>::ToString(ToStringFormat format, const std::string& linesep_in, int digits) const
  {
    char cplx = format == TS_PYTHON ? 'j' : 'i';

    // element separator
    std::string sep = format == TS_PYTHON ? ", " : " ";

    // line separator
    std::string linesep = linesep_in != "" ? linesep_in : (format == TS_PYTHON ? "," : "\n");

    std::ostringstream os;

    // ignored for int case, we cannot fill ints easily as it makes no sense for neg and >= not possible for complex
    if(digits > 0)
      os << std::scientific << std::setprecision(digits);

    switch(format)
    {
    case TS_INFO: {
      int nnz = 0;
      double min = Abs(data_[0][0]);
      double max = Abs(data_[0][0]);
      for(unsigned int k = 0; k < size_row_ * size_col_; k++)
      {
        double v = Abs(data_[0][k]);
        nnz = v != 0 ? nnz+1 : nnz;
        min = std::min(min,v);
        max = std::max(max,v);
      }
      os << "rows=" << size_row_ << " cols=" << size_col_ << " nnz=" << nnz << " min=" << min << " max=" << max;
      break;
    }
    case TS_MATLAB:
    case TS_PYTHON:
      os << "[";
      // intentionally no break;
    case TS_PLAIN:
      for(unsigned int i = 0; i < size_row_; i++)
      {
        if(format == TS_PYTHON)
          os << "[";
        for(unsigned int j = 0; j < size_col_; j++)
        {
          if(boost::is_complex<TYPE>::value)
          {
            Complex cval = (Complex) data_[i][j];
            os << cval.real() << "+" << cval.imag() << cplx;
          }
          else
            os << data_[i][j];

          if(j < size_col_ - 1)
            os << sep;
        }
        if(i < size_row_ - 1)
        {
          if(format == TS_PYTHON)
            os << "]";
          os << linesep;
        }
      }
      if(format == TS_PYTHON)
        os << "]"; // Python closes with ]]
      if(format != TS_PLAIN)
        os << "]";
      break;

    case TS_NONZEROS: {
      for(unsigned int i = 0; i < size_row_; i++)
      {
        bool has_data = false;

        for(unsigned int j = 0; j < size_col_; j++)
        {
          if(Abs(data_[i][j]) != 0.0)
          {
            if(!has_data)
            {
              os << "row=" << i;
              has_data = true;
            }
            os << " " << j << ":";
            if(boost::is_complex<TYPE>::value)
            {
              Complex cval = (Complex) data_[i][j];
              os << cval.real() << "+" << cval.imag() << cplx;
            }
            else
              os << data_[i][j];
          } // end if Abs
        } // end col loop
        if(i < size_row_ - 1 && has_data)
          os << linesep;
      } // end row loop
      break;
    }
    } // end switch

    return os.str();
  }

  /*
    * NOTE: everything (more or less) taken from
    * http://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
    * (3.6.2016)
    * Use upscale to output a 10x10 matrix as an (10*upscale)x(10*upscale) matrix
    *
    * Addition: reference to second matrix which deliveres information for the green channel
    *
    */
    template<class TYPE>
    void Matrix<TYPE>::matrix2Bmp(UInt upscale, std::string filename,Matrix<TYPE>* greenChannel) {
      EXCEPTION("Only implemented for matrices of type double");
    }

    /*
     * Coloring of BMP:
     *  negative values from -1 to 0 are mapped to red color with values from 255 to 0
     *  positive values from 0 to 1 are mapped to blue color with values from 0 to 255
     *
     *  greenChannel > 0: set green channel to 255
     *  greenChallen <= 0 or undefined: set green channel to 0
     */
   template<>
   void Matrix<Double>::matrix2Bmp(UInt upscale, std::string filename,Matrix<Double>* greenChannel)
   {
     if(upscale == 0){
       WARN("Upscaling has to be larger 0");
       return;
     }

     /*
      * Get width and height of matrix and upscale it to image size
      */
     //get dimension of matrix
     UInt height = size_row_*upscale;
     UInt width = size_col_*upscale;

     bool setGreenChannel = false;
     if(greenChannel != NULL){
       if((size_row_ == greenChannel->GetNumRows())&&(size_col_ == greenChannel->GetNumCols())){
         setGreenChannel = true;
       }
     }

     /*
      * Image lines in BMP have to be multiples of 4
      * (*3 because of r g b values)
      */
     UInt padsize = (4-(width*3)%4)%4;
     UInt datasize = (width*3 + padsize) * height;

     /*
      * header and info to be included in BMP files
      */
     unsigned char file[14] = {
         'B','M', // magic
         0,0,0,0, // size in bytes
         0,0, // app data
         0,0, // app data
         40+14,0,0,0 // start of data offset
     };
     unsigned char info[40] = {
         40,0,0,0, // info hd size
         0,0,0,0, // width
         0,0,0,0, // heigth
         1,0, // number color planes
         24,0, // bits per pixel
         0,0,0,0, // compression is none
         0,0,0,0, // image bits size
         0x13,0x0B,0,0, // horz resoluition in pixel / m
         0x13,0x0B,0,0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
         0,0,0,0, // #colors in pallete
         0,0,0,0, // #important colors
         };

     UInt totalsize = datasize + sizeof(file) + sizeof(info);

     /*
      * split all informations into chunks of 1 byte
      */
     file[ 2] = (unsigned char)( totalsize    );
     file[ 3] = (unsigned char)( totalsize>> 8);
     file[ 4] = (unsigned char)( totalsize>>16);
     file[ 5] = (unsigned char)( totalsize>>24);

     info[ 4] = (unsigned char)( width   );
     info[ 5] = (unsigned char)( width>> 8);
     info[ 6] = (unsigned char)( width>>16);
     info[ 7] = (unsigned char)( width>>24);

     info[ 8] = (unsigned char)( height    );
     info[ 9] = (unsigned char)( height>> 8);
     info[10] = (unsigned char)( height>>16);
     info[11] = (unsigned char)( height>>24);

     info[20] = (unsigned char)( datasize    );
     info[21] = (unsigned char)( datasize>> 8);
     info[22] = (unsigned char)( datasize>>16);
     info[23] = (unsigned char)( datasize>>24);

     /*
      * get output stream
      */
     std::ofstream outfile;
     outfile.open(filename.c_str(),std::ofstream::binary);

     if(!outfile.is_open()){
       WARN("Could not open output file!")
       return;
     }

     outfile.write( (char*)file, sizeof(file));
     outfile.write( (char*)info, sizeof(info));

     unsigned char pad[3] = {0,0,0};

     UInt idx, idy;

     /*
      * BMP is stored upside down from right to left
      */
     //for ( UInt y=height-1; y>0; y-- )
     for ( UInt y=0; y<height; y++ )
     {
       //idy = height/upscale-1 - ceil(y/upscale);
       idy = ceil(y/upscale);
       for ( UInt x=0; x<width; x++ )
       {
         //idx = width/upscale-1 - ceil(x/upscale);
         idx = ceil(x/upscale);

         //std::cout << "x, idx: " << x << ", " << idx << std::endl;
         //std::cout << "y, idy: " << y << ", " << idy << std::endl;

         /*
          * Positive values -> blue
          * Negative values -> red
          */
         long red, green, blue;
//         red = lround( -255.0 * data_[idx][idy] );
//         if ( red < 0 ) red=0;
//         if ( red > 255 ) red=255;
//         green = 0;
//         blue = lround( 255.0 * data_[idx][idy] );
//         if ( blue < 0 ) blue=0;
//         if ( blue > 255 ) blue=255;
//         red = 0;
//         green = 0;


         if(data_[idy][idx]<= 0){
           red = lround( -255.0 * data_[idy][idx] );
           if ( red < 0 ) red=0;
           if ( red > 255 ) red=255;
           green = 0;
           blue = 0;
         } else {
           blue = lround( 255.0 * data_[idy][idx] );
           if ( blue < 0 ) blue=0;
           if ( blue > 255 ) blue=255;
           red = 0;
           green = 0;
         }

         if(setGreenChannel){
           /*
            * color only every second pixel (otherwise the figures will have
            * eye-aching colors)
            */
           if((x+y)%2 == 0){
             if((*greenChannel)[idy][idx] > 0){
               /*
                * if we have a positive value of green channel, set green to max of blue and red
                * (otherwise poorly red/blue regions will be only green)
                */
               green = std::max(blue,red);
             } else {
               green = 0;
             }
           } else {
             green = 0;
           }
           //if ( green > 0 ) green=255;
           //if ( green < 0 ) green=0;
         }

         unsigned char pixel[3];
         pixel[0] = blue;
         pixel[1] = green;
         pixel[2] = red;

         outfile.write( (char*)pixel, 3 );
       }
       outfile.write( (char*)pad, padsize );
     }

     outfile.close();
   }

   /*
    * New (29.6.2016)
    *   Outputfunction used for writing out the evaluated state of the vector Preisach model, i.e.
    *   the multiplied rotation and switching state.
    *   Matrix: stores switching state between -1 and 1
    *   RotX: stores x-component of rotation state (between -1 and 1)
    *   RotY: stores y-component of rotation state (between -1 and 1)
    *
    *   Coloring:
    *     - angle between rotation state and x-axis defines color
    *       0 -> red
    *       120 -> blue
    *       240 -> green
    *       0-120 -> red decreases linearly, blue increases linearly, green = 0
    *       120 - 240 -> red = 0, blue decreases linearly, green increases linearly
    *       240 - 360 -> red increases linearly, blue = 0, green decreases linearly
    *     - switching state gives a scaling to the values and a possible offset to the angle
    *       -> switching state < 0 -> 180 degree offset
    *       -> abs(switching state) < 0 -> scale all colors with that value
    * 
    * New concept 11.3.2020
    * - colorrange mapped to 0-180 degree instead of 0-360 degree
    * - for angle between 0-180 degree use:
    *       0 - 60 -> red decreases linearly, blue increases linearly, green = 0
    *       60 - 120 -> red = 0, blue decreases linearly, green increases linearly
    *       120 - 180 -> red increases linearly, blue = 0, green decreases linearly
    * - for angle between -180 and 0 degree set only every second pixel and use:
    *       -180 - -120 -> red decreases linearly, blue increases linearly, green = 0
    *       -120 - -60 -> red = 0, blue decreases linearly, green increases linearly
    *       -60 - 0 -> red increases linearly, blue = 0, green decreases linearly
    * - switching state
    *
   * NOTE: everything (more or less) taken from
   * http://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
   * (3.6.2016)
   * Use upscale to output a 10x10 matrix as an (10*upscale)x(10*upscale) matrix
   *
   * Addition: reference to second matrix which deliveres information for the green channel
   *
   */
   template<class TYPE>
   void Matrix<TYPE>::matrix2Bmp_v2(UInt upscale, std::string filename,Matrix<TYPE>* rotX, Matrix<TYPE>* rotY) {
     EXCEPTION("Only implemented for matrices of type double");
   }

  template<>
  void Matrix<Double>::matrix2Bmp_v2(UInt upscale, std::string filename,Matrix<Double>* rotX, Matrix<Double>* rotY)
  {
    if(upscale == 0){
      WARN("Upscaling has to be larger 0");
      return;
    }

    /*
     * Get width and height of matrix and upscale it to image size
     */
    //get dimension of matrix
    UInt height = size_row_*upscale;
    UInt width = size_col_*upscale;

    if(rotX == NULL || rotY == NULL){
      EXCEPTION("Rotation states are not initialized!");
    }

    /*
     * Image lines in BMP have to be multiples of 4
     * (*3 because of r g b values)
     */
    UInt padsize = (4-(width*3)%4)%4;
    UInt datasize = (width*3 + padsize) * height;

    /*
     * header and info to be included in BMP files
     */
    unsigned char file[14] = {
        'B','M', // magic
        0,0,0,0, // size in bytes
        0,0, // app data
        0,0, // app data
        40+14,0,0,0 // start of data offset
    };
    unsigned char info[40] = {
        40,0,0,0, // info hd size
        0,0,0,0, // width
        0,0,0,0, // heigth
        1,0, // number color planes
        24,0, // bits per pixel
        0,0,0,0, // compression is none
        0,0,0,0, // image bits size
        0x13,0x0B,0,0, // horz resoluition in pixel / m
        0x13,0x0B,0,0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
        0,0,0,0, // #colors in pallete
        0,0,0,0, // #important colors
        };

    UInt totalsize = datasize + sizeof(file) + sizeof(info);

    /*
     * split all informations into chunks of 1 byte
     */
    file[ 2] = (unsigned char)( totalsize    );
    file[ 3] = (unsigned char)( totalsize>> 8);
    file[ 4] = (unsigned char)( totalsize>>16);
    file[ 5] = (unsigned char)( totalsize>>24);

    info[ 4] = (unsigned char)( width   );
    info[ 5] = (unsigned char)( width>> 8);
    info[ 6] = (unsigned char)( width>>16);
    info[ 7] = (unsigned char)( width>>24);

    info[ 8] = (unsigned char)( height    );
    info[ 9] = (unsigned char)( height>> 8);
    info[10] = (unsigned char)( height>>16);
    info[11] = (unsigned char)( height>>24);

    info[20] = (unsigned char)( datasize    );
    info[21] = (unsigned char)( datasize>> 8);
    info[22] = (unsigned char)( datasize>>16);
    info[23] = (unsigned char)( datasize>>24);

    /*
     * get output stream
     */
    std::ofstream outfile;
    outfile.open(filename.c_str(),std::ofstream::binary);

    if(!outfile.is_open()){
      WARN("Could not open output file!")
      return;
    }

    outfile.write( (char*)file, sizeof(file));
    outfile.write( (char*)info, sizeof(info));

    unsigned char pad[3] = {0,0,0};

    UInt idx, idy;
    /*
     * BMP is stored upside down from right to left
     */
    //for ( UInt y=height-1; y>0; y-- )
    for ( UInt y=0; y<height; y++ )
    {
      //idy = height/upscale-1 - ceil(y/upscale);
      idy = ceil(y/upscale);
      for ( UInt x=0; x<width; x++ )
      {
        //idx = width/upscale-1 - ceil(x/upscale);
        idx = ceil(x/upscale);

        long red, green, blue;

        Double scaling = abs(data_[idy][idx]);

        /*
         * OLD
         *  On diagonal idx = idy we have half value
         *    -> color magnitude halved
         *
         * NEW (22.8.16)
         *  On diagonal idx = idy, double value but color only the pixels y>=x
         *
         */

        if(idx == idy){
          if(y >= x){
            scaling = 2.0*scaling;
          }
//          else {
//            /*
//             * make pixel white
//             */
//            red = 255;
//            blue = 255;
//            green = 255;
//          }
        }

        /*
         * calculate angle for determination of coloring
         */
        Double tmp = (*rotX)[idy][idx];
        if(tmp > 1.0){
          tmp = 1.0;
        } else if(tmp < -1.0){
          tmp = -1.0;
        }
        Double angle = std::acos(tmp)/M_PI * 180;
        if((*rotY)[idy][idx] < 0){
          angle = 360-angle;
        }
        if(data_[idy][idx] < 0){
          angle = 180+angle;
        }
        if(angle >= 360){
          angle = angle - 360;
        }

//        if(x == 1 && y == height-3){
//          std::cout << "x,y: " << x << ", " << y << std::endl;
//          std::cout << "rotX,rotY: " << (*rotX)[idy][idx] << ", " << (*rotY)[idy][idx] << std::endl;
//          std::cout << "scaling,angle: " << scaling << ", " << angle << std::endl;
//        }
//
//        if(x == width-3 && y == height-3){
//          std::cout << "x,y: " << x << ", " << y << std::endl;
//          std::cout << "rotX,rotY: " << (*rotX)[idy][idx] << ", " << (*rotY)[idy][idx] << std::endl;
//          std::cout << "scaling,angle: " << scaling << ", " << angle << std::endl;
//        }
//

        if(angle >= 0 && angle < 60){
          red = lround(scaling*255.0);
        } else if(angle >= 60 && angle < 120){
          red = lround(scaling*255.0*(120-angle)/60.0);
        } else if(angle >= 240 && angle < 300){
          red = lround(scaling*255.0*(angle-240)/60.0);
        } else if(angle >= 300 && angle < 360){
          red = lround(scaling*255.0);
        } else {
          red = 0;
        }

        if(angle >= 60 && angle < 180){
          blue = lround(scaling*255.0);
        } else if(angle >= 180 && angle < 240){
          blue = lround(scaling*255.0*(240-angle)/60.0);
        } else if(angle >= 0 && angle < 60){
          blue = lround(scaling*255.0*angle/60.0);
        } else {
          blue = 0;
        }

        if(angle >= 180 && angle < 300){
          green = lround(scaling*255.0);
        } else if(angle >= 120 && angle < 180){
          green = lround(scaling*255.0*(angle-120)/60.0);
        } else if(angle >= 300 && angle < 360){
          green = lround(scaling*255.0*(360-angle)/60.0);
        } else {
          green = 0;
        }


//        if(x == lround(width/4.0) && y == lround(0.6*height)){
//          std::cout << "x,y: " << x << ", " << y << std::endl;
//          std::cout << "rotX,rotY: " << (*rotX)[idy][idx] << ", " << (*rotY)[idy][idx] << std::endl;
//          std::cout << "scaling,angle: " << scaling << ", " << angle << std::endl;
//          std::cout << "rgb: " << red << "," << green << "," << blue << std::endl;
//        }
//        if(x == 1){
//          std::cout << "rotX, rotY: " << (*rotX)[idy][idx] << "," << (*rotY)[idy][idx] << std::endl;
//          std::cout << "rgb: " << red << "," << green << "," << blue << std::endl;
//        }

        if((*rotX)[idy][idx] == 0 && (*rotY)[idy][idx] == 0){
//          /*
//           * no rotation state -> make every second pixel white!
//           */
//          if((x+y)%2 == 0){
//            red = 255;
//            blue = 255;
//            green = 255;
//          }
          /*
           * no rotation state -> make every pixel white, but only for point below the diagonal!
           */
          if(y >= x){
            /*
             * -> above alpha = beta
             */
            // every pixel gray
            red = 100;
            blue = 100;
            green = 100;

            if(x+y+1 > height){
              /*
               * above diagonal alpha = -beta
               */
              // every second pixel black
              if((x+y)%2 == 0){
                red = 0;
                blue = 0;
                green = 0;
              }
            }
          } else {
            /*
             * every pixel white
             */
            red = 255;
            blue = 255;
            green = 255;
          }
        }

//
//        if(angle >= 0 && angle < 120){
//          red = lround(scaling*255.0*(120-angle)/120.0);
//          blue = lround(scaling*255.0*angle/120.0);
//          green = 0;
//        } else if(angle >= 120 && angle < 240){
//          red = 0;
//          blue = lround(scaling*255.0*(240-angle)/120.0);
//          green = lround(scaling*255.0*(angle-120)/120.0);
//        } else if(angle >= 240 && angle < 360){
//          red = lround(scaling*255.0*(360-angle)/120.0);
//          blue = 0;
//          green = lround(scaling*255.0*(angle-240)/120.0);l
//        } else {
//          WARN("This angle should not occur!");
//        }

        if ((*rotX)[idy][idx] != 0 || (*rotY)[idy][idx] != 0){
          if(data_[idy][idx] < 0){
            /*
             * new coloring: angles > 180 have same color as angle-180 but only every second pixel is colored
             */
            // every second pixel white
            if((x+y)%2 == 0){
              red = 0;
              blue = 0;
              green = 0;
            } else {
              // every other pixel flips color
              red = 255-red;
              blue = 255-blue;
              green = 255-green;
            }
          }
        }

        if(idx == idy){
          if(y < x){
            /*
             * make pixel white in lower half of main diagonal to form triangles
             */
            red = 255;
            blue = 255;
            green = 255;
          }
        }

        unsigned char pixel[3];
        pixel[0] = blue;
        pixel[1] = green;
        pixel[2] = red;

        outfile.write( (char*)pixel, 3 );
      }
      outfile.write( (char*)pad, padsize );
    }

    outfile.close();
  } 
   
  void GetColor_from_ColorMap(UInt colormapId, Double value, Double brightness, long* red, long* green, long* blue){

    if(colormapId == 1){
      /*
       * diverging colormap
       */
      *red = lround( brightness*(value/360.0));
      *green = 0;
      *blue = lround( brightness*((360-value)/360.0));
      
    } else {
      /*
       * 360 degree cyclic colormap
       */
      if(value >= 0 && value < 120){
        *red = lround( brightness*(120-value)/120.0 ) ;
        *blue = lround( brightness*(value/120.0 ) );
        *green = 0;
      } else if(value >= 120 && value < 240){
        *red = 0;
        *blue = lround( brightness*(240-value)/120.0 );
        *green = lround( brightness*(value-120)/120.0 );
      } else {
        *red = lround( brightness*(value-240)/120.0 );
        *blue = 0;
        *green = lround( brightness*(360-value)/120.0 );
      }
    }
  }
  
  template<class TYPE>
   void Matrix<TYPE>::matrix2Bmp_v3(UInt upscale, std::string filename,Matrix<TYPE>* rotX, Matrix<TYPE>* rotY) {
     EXCEPTION("Only implemented for matrices of type double");
   }
  
  template<>
  void Matrix<Double>::matrix2Bmp_v3(UInt upscale, std::string filename,Matrix<Double>* rotX, Matrix<Double>* rotY) {
    
    if(upscale == 0){
      WARN("Upscaling has to be larger 0");
      return;
    }

    /*
     * Get width and height of matrix and upscale it to image size
     */
    //get dimension of matrix
    UInt height = size_row_*upscale;
    UInt width = size_col_*upscale;
    
    /*
     * Extension 13.3.2020
     * - Previously rotX and rotY had to be specified, i.e., this function crashed if no rotation state was delivered
     *   In consequence, this function only was usable to output an overlay of a 2d rotation state with the data of 
     *   this matrix. 
     * - To make this function more useful, e.g., to print out only the data which could represent the switching state
     *   of the Preisach plane or the weight distribution function, the following switch is introduced:
     *     rotStateGiven
     * - If rotStateGiven = false, the coloring of cell i,j will no longer be according to the angle atan2(rotY[i][j],rotX[i][j])
     *   but according to the normalized value of data[i][j]; this leads also to different annotations and labels.
     *   Note however, that still an upper triangular shape is assuemed!
     * 
     */
    bool rotStateGiven = true;
    UInt colormapId = 2; //2 = cyclic colormap
    Double data_abs_max = 1.0;
    
    if(rotX == NULL || rotY == NULL){
//      EXCEPTION("Rotation states are not initialized!");
      rotX = this;
      rotY = this;
      rotStateGiven = false;
      data_abs_max = std::max(this->GetMax(),-this->GetMin());
      colormapId = 1; // 1 = diverging colormap
    }

    /*
     * Image lines in BMP have to be multiples of 4
     * (*3 because of r g b values)
     */
    UInt padsize = (4-(width*3)%4)%4;
    UInt datasize = (width*3 + padsize) * height;

    /*
     * header and info to be included in BMP files
     */
    unsigned char file[14] = {
        'B','M', // magic
        0,0,0,0, // size in bytes
        0,0, // app data
        0,0, // app data
        40+14,0,0,0 // start of data offset
    };
    unsigned char info[40] = {
        40,0,0,0, // info hd size
        0,0,0,0, // width
        0,0,0,0, // heigth
        1,0, // number color planes
        24,0, // bits per pixel
        0,0,0,0, // compression is none
        0,0,0,0, // image bits size
        0x13,0x0B,0,0, // horz resoluition in pixel / m
        0x13,0x0B,0,0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
        0,0,0,0, // #colors in pallete
        0,0,0,0, // #important colors
        };

    UInt totalsize = datasize + sizeof(file) + sizeof(info);

    /*
     * split all informations into chunks of 1 byte
     */
    file[ 2] = (unsigned char)( totalsize    );
    file[ 3] = (unsigned char)( totalsize>> 8);
    file[ 4] = (unsigned char)( totalsize>>16);
    file[ 5] = (unsigned char)( totalsize>>24);

    info[ 4] = (unsigned char)( width   );
    info[ 5] = (unsigned char)( width>> 8);
    info[ 6] = (unsigned char)( width>>16);
    info[ 7] = (unsigned char)( width>>24);

    info[ 8] = (unsigned char)( height    );
    info[ 9] = (unsigned char)( height>> 8);
    info[10] = (unsigned char)( height>>16);
    info[11] = (unsigned char)( height>>24);

    info[20] = (unsigned char)( datasize    );
    info[21] = (unsigned char)( datasize>> 8);
    info[22] = (unsigned char)( datasize>>16);
    info[23] = (unsigned char)( datasize>>24);

    /*
     * get output stream
     */
    std::ofstream outfile;
    outfile.open(filename.c_str(),std::ofstream::binary);

    if(!outfile.is_open()){
      WARN("Could not open output file!")
      return;
    }

    outfile.write( (char*)file, sizeof(file));
    outfile.write( (char*)info, sizeof(info));

    unsigned char pad[3] = {0,0,0};

    UInt idx, idy;
    
    /*
     * Colorbar for rotation state
     */
    UInt colorBarWidth = ceil((0.05*width)/upscale);
    UInt colorBarHeight = ceil((0.5*height)/upscale);
    UInt colorBarXOffset = ceil((0.8*width)/upscale);
    UInt colorBarYOffset = ceil((0.1*height)/upscale);
    UInt numColors = 360;
    UInt colorBarHeightUnsetState = ceil(1.5*colorBarWidth);
//    std::cout << "Colorbar-Info:" << std::endl;
//    std::cout << "colorBarWidth in px: " << colorBarWidth << std::endl;
//    std::cout << "colorBarHeight in px: " << colorBarHeight << std::endl;
//    std::cout << "colorBarXOffset in px: " << colorBarXOffset << std::endl;
//    std::cout << "colorBarYOffset in px: " << colorBarYOffset << std::endl;
//    std::cout << "numberOfColors: " << numColors << std::endl;
//    
    /*
     * Define crude bitmaps for labels
     */
    bool triggerLabel = false;
    UInt labelCnt = 0;
    UInt labelXOffset = colorBarXOffset + colorBarWidth;
    UInt labelYOffsetBase = colorBarYOffset;
    UInt labelYOffset = labelYOffsetBase;
    UInt deltaYOffset = round(colorBarHeight/4.0);
    
    static UInt labelBaseWidth = 24;
    static UInt labelBaseHeight = 8;
    
    Double upscaleLabel = 3*upscale;
    UInt labelWidth = labelBaseWidth*upscaleLabel;
    UInt labelHeight = labelBaseHeight*upscaleLabel;
    UInt titleAddYOffset = ceil(labelHeight/2.0);
    UInt titleXOffset = colorBarXOffset;

//    std::cout << "upscaleLabel = " << upscaleLabel << std::endl;
//    std::cout << "labelWidth = " << labelWidth << std::endl;
//    std::cout << "labelHeight = " << labelHeight << std::endl;

    std::vector<int>* bmpNeg180 = new std::vector<int>[labelBaseHeight];
    std::vector<int>* bmpNeg90 = new std::vector<int>[labelBaseHeight];
    std::vector<int>* bmp0 = new std::vector<int>[labelBaseHeight];
    std::vector<int>* bmpPos90 = new std::vector<int>[labelBaseHeight];
    std::vector<int>* bmpPos180 = new std::vector<int>[labelBaseHeight];
    std::vector<int>* colorbarTitle = new std::vector<int>[labelBaseHeight];
    std::vector<int>* unsetRotState = new std::vector<int>[labelBaseHeight];
    
    if(rotStateGiven){
      /*
       * colorbar corresponds to angle; bmp maps / labels shall indicate angle
       */
      bmpNeg180[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpNeg180[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpNeg180[2] = {0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpNeg180[3] = {0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0};
      bmpNeg180[4] = {0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,0,1,0,1,0,0,0,0};
      bmpNeg180[5] = {0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0};
      bmpNeg180[6] = {0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpNeg180[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmpNeg90[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpNeg90[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpNeg90[2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpNeg90[3] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0};
      bmpNeg90[4] = {0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,1,0,1,0,0,0,0};
      bmpNeg90[5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0};
      bmpNeg90[6] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpNeg90[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmp0[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmp0[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmp0[2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0};
      bmp0[3] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0};
      bmp0[4] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0};
      bmp0[5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0};
      bmp0[6] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0};
      bmp0[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmpPos90[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpPos90[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpPos90[2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpPos90[3] = {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0};
      bmpPos90[4] = {0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,1,0,1,0,0,0,0};
      bmpPos90[5] = {0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0};
      bmpPos90[6] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpPos90[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmpPos180[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpPos180[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpPos180[2] = {0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpPos180[3] = {0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0};
      bmpPos180[4] = {0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,0,1,0,1,0,0,0,0};
      bmpPos180[5] = {0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0};
      bmpPos180[6] = {0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,0,0,0};
      bmpPos180[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      colorbarTitle[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      colorbarTitle[1] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      colorbarTitle[2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      colorbarTitle[3] = {0,1,0,1,0,1,1,1,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0};
      colorbarTitle[4] = {0,1,0,1,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,1,0};
      colorbarTitle[5] = {0,1,1,1,0,1,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0};
      colorbarTitle[6] = {0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0};
      colorbarTitle[7] = {0,1,1,1,0,1,0,0,0,1,1,1,0,1,0,1,0,1,1,1,0,0,0,0};
      
      unsetRotState[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      unsetRotState[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      unsetRotState[2] = {0,0,0,1,1,1,0,1,0,0,1,0,1,1,1,0,1,1,1,0,0,1,0,0};
      unsetRotState[3] = {0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0};
      unsetRotState[4] = {0,0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0};
      unsetRotState[5] = {0,0,0,1,0,1,0,1,1,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0};
      unsetRotState[6] = {0,0,0,1,0,1,0,1,0,0,1,0,1,1,1,0,1,1,1,0,1,1,1,0};
      unsetRotState[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
    } else {
      /*
       * No rotstate given; colorbar shall represent normalized amplitude of data field
       * keep name of bmp to reuse algorithm but change pattern of 
       */
      /*
       * -180 degree -> amplitdue -1,0
       * -90 degree -> amplitude -0.5 
       * ...
       */
      bmpNeg180[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpNeg180[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpNeg180[2] = {0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,1,1,0,0,0};
      bmpNeg180[3] = {0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmpNeg180[4] = {0,0,0,0,0,1,1,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmpNeg180[5] = {0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmpNeg180[6] = {0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,0};
      bmpNeg180[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmpNeg90[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpNeg90[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpNeg90[2] = {0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,0};
      bmpNeg90[3] = {0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,0};
      bmpNeg90[4] = {0,0,0,0,0,1,1,1,0,1,0,1,0,0,1,1,1,0,1,0,1,0,0,0};
      bmpNeg90[5] = {0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0};
      bmpNeg90[6] = {0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0};
      bmpNeg90[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmp0[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmp0[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmp0[2] = {0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,0};
      bmp0[3] = {0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmp0[4] = {0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmp0[5] = {0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmp0[6] = {0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0};
      bmp0[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmpPos90[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpPos90[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpPos90[2] = {0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,0};
      bmpPos90[3] = {0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,0};
      bmpPos90[4] = {0,0,0,0,0,1,1,1,0,1,0,1,0,0,1,1,1,0,1,0,1,0,0,0};
      bmpPos90[5] = {0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0};
      bmpPos90[6] = {0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0};
      bmpPos90[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      bmpPos180[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      bmpPos180[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      bmpPos180[2] = {0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,1,1,0,0,0};
      bmpPos180[3] = {0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmpPos180[4] = {0,0,0,0,0,1,1,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmpPos180[5] = {0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0};
      bmpPos180[6] = {0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,0};
      bmpPos180[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
      /*
       * AMPL. instead of ALPHA
       */
      colorbarTitle[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      colorbarTitle[1] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      colorbarTitle[2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      colorbarTitle[3] = {0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,0};
      colorbarTitle[4] = {0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0};
      colorbarTitle[5] = {0,1,1,1,0,1,0,1,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0};
      colorbarTitle[6] = {0,1,0,1,0,1,1,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,1,0};
      colorbarTitle[7] = {0,1,1,1,0,1,0,0,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0};
      
      unsetRotState[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      unsetRotState[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      unsetRotState[2] = {0,0,0,1,1,1,0,1,0,0,1,0,1,1,1,0,1,1,1,0,0,1,0,0};
      unsetRotState[3] = {0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0};
      unsetRotState[4] = {0,0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0};
      unsetRotState[5] = {0,0,0,1,0,1,0,1,1,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0};
      unsetRotState[6] = {0,0,0,1,0,1,0,1,0,0,1,0,1,1,1,0,1,1,1,0,1,1,1,0};
      unsetRotState[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    }
    std::vector<int>* currentLabel = bmpNeg180;
    
    /*
     * Additional labels for switching states
     * - bright color -> +1
     * - shaded color -> -1
     */
    bool triggerLabel2 = false;
    UInt labelCnt2 = 0;
    UInt colorBarXOffset2 = ceil(0.7*colorBarXOffset);
    UInt colorBarYOffset2 = colorBarYOffset;
    UInt colorBarWidth2 = 3*colorBarWidth;
    UInt labelXOffset2 = colorBarXOffset2 + colorBarWidth2;
    UInt labelYOffset2 = colorBarYOffset2 + ceil(0.5*colorBarWidth);
    
    static UInt labelBaseWidth2 = 29;
    static UInt labelBaseHeight2 = 8;
    
    UInt labelWidth2 = labelBaseWidth2*upscaleLabel;
    UInt labelHeight2 = labelBaseHeight2*upscaleLabel;
    UInt titleAddYOffset2 = ceil(labelHeight/2.0);

    std::vector<int>* switchingTitle = new std::vector<int>[labelBaseHeight2];
    switchingTitle[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    switchingTitle[1] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    switchingTitle[2] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    switchingTitle[3] = {0,1,1,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,0,1,1,1,0,1,0,1,0,0,0};
    switchingTitle[4] = {0,0,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,1,0};
    switchingTitle[5] = {0,1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,1,0,0,1,0,0,0,1,1,1,0,0,0};
    switchingTitle[6] = {0,1,0,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,1,0};
    switchingTitle[7] = {0,1,1,1,0,1,0,0,0,0,0,1,0,1,0,1,1,1,0,1,1,1,0,1,0,1,0,0,0};
        
    std::vector<int>* plus1 = new std::vector<int>[labelBaseHeight2];
    plus1[0] = {1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    plus1[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    plus1[2] = {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    plus1[3] = {0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    plus1[4] = {0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    plus1[5] = {0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    plus1[6] = {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    plus1[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    std::vector<int>* minus1 = new std::vector<int>[labelBaseHeight2];
    minus1[0] = {1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    minus1[1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    minus1[2] = {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    minus1[3] = {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    minus1[4] = {0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    minus1[5] = {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    minus1[6] = {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    minus1[7] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    std::vector<int>* currentLabel2 = switchingTitle;
    
    /*
     * for grid
     * > plot lines at -1.0 / -0.8 / -0.6 / -0.4 / -0.2 / 0.0 / 0.2 / 0.4 / 0.6 / 0.8 / 1.0
     */
    bool printGrid = true;
    bool triggerLabel3 = false;
    UInt labelCnt3 = 0;
    static UInt labelBaseWidth3 = 16;
    static UInt labelBaseHeight3 = 16;
    
    UInt labelWidth3 = labelBaseWidth3*upscaleLabel;
    UInt labelHeight3 = labelBaseHeight3*upscaleLabel;
    UInt blockWidth = labelWidth3;
    UInt blockHeight = blockWidth;        
    UInt blockXOffset = ceil(0.4*colorBarXOffset);
    UInt blockYOffset = colorBarYOffset;
    UInt labelXOffset3 = blockXOffset + blockWidth;
    UInt labelYOffset3 = blockYOffset + blockHeight;
        
    // for ten cells we have a cell size of 0.2 x 0.2 
    // > if number of lines should be changed define a fitting bitmap
    std::vector<int>* zeroDotTwoTop = new std::vector<int>[labelBaseHeight3];
    zeroDotTwoTop[0] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    zeroDotTwoTop[1] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    zeroDotTwoTop[2] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    zeroDotTwoTop[3] = {0,0,0,0,1,1,1,1,0,1,1,1,0,0,0,0};
    zeroDotTwoTop[4] = {0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0};
    zeroDotTwoTop[5] = {0,0,0,0,1,0,1,0,0,1,1,1,0,0,0,0};
    zeroDotTwoTop[6] = {0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0};
    zeroDotTwoTop[7] = {0,0,0,0,1,1,1,0,0,1,1,1,0,0,0,0};
    zeroDotTwoTop[8] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoTop[9] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoTop[10] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoTop[11] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoTop[12] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoTop[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoTop[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoTop[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    std::vector<int>* zeroDotTwoSide = new std::vector<int>[labelBaseHeight3];
    zeroDotTwoSide[0] = {1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[1] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[2] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[3] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[4] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[5] = {0,1,0,1,1,1,1,0,1,1,1,0,0,0,0,0};
    zeroDotTwoSide[6] = {0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0};
    zeroDotTwoSide[7] = {0,1,0,1,0,1,0,0,1,1,1,0,0,0,0,0};
    zeroDotTwoSide[8] = {0,1,0,1,0,1,0,0,0,0,1,0,0,0,0,0};
    zeroDotTwoSide[9] = {0,1,0,1,1,1,0,0,1,1,1,0,0,0,0,0};
    zeroDotTwoSide[10] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[11] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[12] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[13] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[14] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    zeroDotTwoSide[15] = {1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    std::vector<int>* currentLabel3 = zeroDotTwoTop;
    
    /*
     * BMP is stored upside down from right to left
     */
    //for ( UInt y=height-1; y>0; y-- )
    long red, green, blue;
    Double scaling;
    
    for ( UInt y=0; y<height; y++ ) {
      
      //idy = height/upscale-1 - ceil(y/upscale);
      idy = ceil(y/upscale);
      for ( UInt x=0; x<width; x++ ) {
        
        //idx = width/upscale-1 - ceil(x/upscale);
        idx = ceil(x/upscale);

        scaling = abs(data_[idy][idx]);

        /*
         *  On diagonal idx = idy, double value but color only the pixels y>=x
         */

        if(idx == idy){
          if(y >= x){
            scaling = 2.0*scaling;
          }
        }

        /*
         * calculate angle for determination of coloring
         */       
        Double rot_x = (*rotX)[idy][idx];
        Double rot_y = (*rotY)[idy][idx];
        
        if(rot_x == 0 && rot_y == 0){
          /*
           * no rotation state 
           * -> below main diagonal: make every pixel white (to get triangular Preisach-plane)
           * -> at right hand side, add colorbar
           * -> above main diagonal: grey if switching state is +1; grey/black if switching state is -1
           */
          if(y >= x){
            /*
             * -> above alpha = beta
             */
            // every pixel gray
            red = 100;
            blue = 100;
            green = 100;

            if(x+y+1 > height){
              /*
               * above diagonal alpha = -beta 
               * > initial switching state = -1
               */
              // every second pixel black
              if((x+y)%2 == 0){
                red = 0;
                blue = 0;
                green = 0;
              }
            }
          } else {
            /*
             * every pixel white
             */
            red = 255;
            blue = 255;
            green = 255;
            
            /*
             * print colorbar
             */
            Double deltaAnglePerColor = 360.0/( (Double) numColors );
            Double pixelPerColor = ((Double) colorBarHeight)/((Double) numColors);
            Double currentAngle = 0.0;
            UInt colorIndex = 0;
            Double magnitude = 255.0;
            
            if((colorBarYOffset <= y) && (y < colorBarYOffset+colorBarHeight)){
              if((colorBarXOffset <= x) && (x < colorBarXOffset+colorBarWidth)){

                colorIndex = floor( (y-colorBarYOffset)/(pixelPerColor));
                // revert colorbar such that 360 degree is on top and 0 degree on bottom
                // remember: - angle is coloring of Preisach plane is shifted by 180 degree, i.e., 
                //            alignment along pos. x-axis results in 180
                //            alignment along pos. y-axis results in 270
                //           - colorbar should indicate the correct angle though, i.e., when adding label
                //            we have 0 degree in the middle and +/-180 degree on the ends
                currentAngle = deltaAnglePerColor*colorIndex;
                
                GetColor_from_ColorMap(colormapId, currentAngle, magnitude, &red, &green, &blue); 
                
//                if(currentAngle >= 0 && currentAngle < 120){
//                  red = lround( magnitude*(120-currentAngle)/120.0 ) ;
//                  blue = lround( magnitude*(currentAngle/120.0 ) );
//                  green = 0;
//                } else if(currentAngle >= 120 && currentAngle < 240){
//                  red = 0;
//                  blue = lround( magnitude*(240-currentAngle)/120.0 );
//                  green = lround( magnitude*(currentAngle-120)/120.0 );
//                } else {
//                  red = lround( magnitude*(currentAngle-240)/120.0 );
//                  blue = 0;
//                  green = lround( magnitude*(360-currentAngle)/120.0 );
//                }

                /*
                 * trigger labels
                 */
                if(((currentAngle <= 0)&&(currentAngle > 0-deltaAnglePerColor))&&(triggerLabel==false)){
//                  std::cout << "label 1 triggered" << std::endl;
                  triggerLabel = true;
                  labelCnt = 0;
                  currentLabel = bmpNeg180;
                  labelYOffset = labelYOffsetBase;
//                  std::cout << "labelYOffset = " << labelYOffset << std::endl;
                }
                if(((currentAngle <= 90)&&(currentAngle > 90-deltaAnglePerColor))&&(triggerLabel==false)){
//                  std::cout << "label 2 triggered" << std::endl;
                  triggerLabel = true;
                  labelCnt = 0;
                  currentLabel = bmpNeg90;
                  labelYOffset = floor(labelYOffset+deltaYOffset);
//                  std::cout << "labelYOffset = " << labelYOffset << std::endl;
                }
                if(((currentAngle <= 180)&&(currentAngle > 180-deltaAnglePerColor))&&(triggerLabel==false)){
//                  std::cout << "label 3 triggered" << std::endl;
                  triggerLabel = true;
                  labelCnt = 0;
                  currentLabel = bmp0;
                  labelYOffset = floor(labelYOffset+deltaYOffset);
//                  std::cout << "labelYOffset = " << labelYOffset << std::endl;
                }
                if(((currentAngle < 270)&&(currentAngle >= 270-deltaAnglePerColor))&&(triggerLabel==false)){
//                  std::cout << "label 4 triggered" << std::endl;
                  triggerLabel = true;
                  labelCnt = 0;
                  currentLabel = bmpPos90;
                  labelYOffset = floor(labelYOffset+deltaYOffset);
//                  std::cout << "labelYOffset = " << labelYOffset << std::endl;
                }
                if(((currentAngle < 360)&&(currentAngle >= 360-deltaAnglePerColor))&&(triggerLabel==false)){
//                  std::cout << "label 5 triggered" << std::endl;
                  triggerLabel = true;
                  labelCnt = 0;
                  currentLabel = bmpPos180;
                  labelYOffset = floor(labelYOffset+deltaYOffset);
//                  std::cout << "labelYOffset = " << labelYOffset << std::endl;
                }
              }
            }
            
            if((ceil(colorBarYOffset+colorBarHeight+1.0/3.0*colorBarHeightUnsetState) <= y) && (y < colorBarYOffset+colorBarHeight+colorBarHeightUnsetState)){
              if((colorBarXOffset <= x) && (x < colorBarXOffset+colorBarWidth)){
                /*
                 * grey block for unset state
                 */
                red = 100;
                green = 100;
                blue = 100;
                
                /*
                 * trigger labels
                 */
                if(y <= (ceil(colorBarYOffset+colorBarHeight+1.0/3.0*colorBarHeightUnsetState) + 1) && (triggerLabel==false)){
                  triggerLabel = true;
                  labelCnt = 0;
                  currentLabel = unsetRotState;
                  labelYOffset = ceil(colorBarYOffset+colorBarHeight+2.0/3.0*colorBarHeightUnsetState);
                }
                
              }
            }
                        
            if((colorBarYOffset+colorBarHeight+colorBarHeightUnsetState+titleAddYOffset <= y) && (y < colorBarYOffset+colorBarHeight+colorBarHeightUnsetState+2*titleAddYOffset)){
              if((colorBarXOffset <= x) && (x < colorBarXOffset+colorBarWidth)){
                if((y <= colorBarYOffset+colorBarHeight+colorBarHeightUnsetState+titleAddYOffset + 1) && (triggerLabel==false)){
                  triggerLabel = true;
                  labelCnt = 0;
                  currentLabel = colorbarTitle;
                  labelYOffset = colorBarYOffset+colorBarHeight+colorBarHeightUnsetState+titleAddYOffset;
                  labelXOffset = titleXOffset;
                } 
              }
            }
            
            /*
             * for switching state > this label is only required in addition to the rotation state
             * > if colorbar indicates the amplitude, simply skip this additional labels
             */
            if(rotStateGiven){
              if((colorBarYOffset2 <= y) && (y < colorBarYOffset2+colorBarWidth)){
                if((colorBarXOffset2 <= x) && (x < colorBarXOffset2+colorBarWidth2)){

                  if((colorBarXOffset2 <= x) && (x < colorBarXOffset2+ceil(colorBarWidth2/5.0))){
                    currentAngle = 0;                  
                  }
                  else if((colorBarXOffset2+ceil(colorBarWidth2/5.0) <= x) && (x < colorBarXOffset2+ceil(2*colorBarWidth2/5.0))){
                    currentAngle = 90;                  
                  }
                  else if((colorBarXOffset2+ceil(2*colorBarWidth2/5.0) <= x) && (x < colorBarXOffset2+ceil(3*colorBarWidth2/5.0))){
                    currentAngle = 180;                  
                  }
                  else if((colorBarXOffset2+ceil(3*colorBarWidth2/5.0) <= x) && (x < colorBarXOffset2+ceil(4*colorBarWidth2/5.0))){
                    currentAngle = 270;                  
                  }
                  else {
                    currentAngle = 360;                  
                  }

                  GetColor_from_ColorMap(colormapId, currentAngle, magnitude, &red, &green, &blue); 
//                  
//                  if(divergingColorMap){
//                    /*
//                     * diverging colormap
//                     */
//                    red = lround( magnitude*(currentAngle/360.0));
//                    green = 0;
//                    blue = lround( magnitude*((360-currentAngle)/360.0));
//                    
//                  } else {
//                    /*
//                     * 360 degree cyclic colormap
//                     */
//                    if(currentAngle >= 0 && currentAngle < 120){
//                      red = lround( magnitude*(120-currentAngle)/120.0 ) ;
//                      blue = lround( magnitude*(currentAngle/120.0 ) );
//                      green = 0;
//                    } else if(currentAngle >= 120 && currentAngle < 240){
//                      red = 0;
//                      blue = lround( magnitude*(240-currentAngle)/120.0 );
//                      green = lround( magnitude*(currentAngle-120)/120.0 );
//                    } else {
//                      red = lround( magnitude*(currentAngle-240)/120.0 );
//                      blue = 0;
//                      green = lround( magnitude*(360-currentAngle)/120.0 );
//                    }
//                  }


                  // every second pixel black
                  if((x+y)%2 == 0){
                    red = 0;
                    blue = 0;
                    green = 0;
                  } 

                  if((y <= colorBarYOffset2 + 1) && (triggerLabel2==false)){
                    triggerLabel2 = true;
                    labelCnt2 = 0;
                    currentLabel2 = minus1;
  //                  labelYOffset2 = labelYOffset2;
  //                  labelXOffset2 = labelXOffset2;
                  } 
                }
              }

              if(((colorBarYOffset2+ceil(1.5*colorBarWidth) <= y)) && (y < colorBarYOffset2+ceil(2.5*colorBarWidth))){
                if((colorBarXOffset2 <= x) && (x < colorBarXOffset2+colorBarWidth2)){

                  if((colorBarXOffset2 <= x) && (x < colorBarXOffset2+ceil(colorBarWidth2/5.0))){
                    currentAngle = 0;                  
                  }
                  else if((colorBarXOffset2+ceil(colorBarWidth2/5.0) <= x) && (x < colorBarXOffset2+ceil(2*colorBarWidth2/5.0))){
                    currentAngle = 90;                  
                  }
                  else if((colorBarXOffset2+ceil(2*colorBarWidth2/5.0) <= x) && (x < colorBarXOffset2+ceil(3*colorBarWidth2/5.0))){
                    currentAngle = 180;                  
                  }
                  else if((colorBarXOffset2+ceil(3*colorBarWidth2/5.0) <= x) && (x < colorBarXOffset2+ceil(4*colorBarWidth2/5.0))){
                    currentAngle = 270;                  
                  }
                  else {
                    currentAngle = 360;                  
                  }

                  GetColor_from_ColorMap(colormapId, currentAngle, magnitude, &red, &green, &blue); 
//                  
//                  if(currentAngle >= 0 && currentAngle < 120){
//                    red = lround( magnitude*(120-currentAngle)/120.0 ) ;
//                    blue = lround( magnitude*(currentAngle/120.0 ) );
//                    green = 0;
//                  } else if(currentAngle >= 120 && currentAngle < 240){
//                    red = 0;
//                    blue = lround( magnitude*(240-currentAngle)/120.0 );
//                    green = lround( magnitude*(currentAngle-120)/120.0 );
//                  } else {
//                    red = lround( magnitude*(currentAngle-240)/120.0 );
//                    blue = 0;
//                    green = lround( magnitude*(360-currentAngle)/120.0 );
//                  }

                  if((y <= colorBarYOffset2+ceil(1.5*colorBarWidth) + 1) && (triggerLabel2==false)){
                    triggerLabel2 = true;
                    labelCnt2 = 0;
                    currentLabel2 = plus1;
                    labelYOffset2 = labelYOffset2+ceil(1.5*colorBarWidth);
  //                  labelXOffset2 = labelXOffset2;
                  } 
                }
              }

              if((colorBarYOffset2+ceil(2.5*colorBarWidth)+titleAddYOffset2 <= y) 
                      && (y < colorBarYOffset2+ceil(2.5*colorBarWidth)+2*titleAddYOffset2)){
                if((colorBarXOffset2 <= x) && (x < colorBarXOffset2+colorBarWidth2)){

                  if((y <= colorBarYOffset2+ceil(2.5*colorBarWidth)+titleAddYOffset2 + 1) && (triggerLabel2==false)){
                    triggerLabel2 = true;
                    labelCnt2 = 0;
                    currentLabel2 = switchingTitle;
                    labelYOffset2 = colorBarYOffset2+ceil(2.5*colorBarWidth)+titleAddYOffset2;
                    labelXOffset2 = colorBarXOffset2;
                  } 
                }
              }
            }
            
            /*
             * write out labels
             */
            if(triggerLabel == true){
              
              if((labelYOffset <= y) && (y < labelYOffset+labelHeight)){
                UInt idRow = floor((y-labelYOffset)/upscaleLabel);
              
                if((labelXOffset <= x) && (x < labelXOffset+labelWidth)){
                  UInt idCol = floor((x-labelXOffset)/upscaleLabel);
//                  std::cout << "idCol, idRow: " << idCol << ", " << idRow << std::endl;
                  
                  // go over bitmap; 1 = black pixel
                  if(currentLabel[idRow][idCol] == 1){
                    red = 0;
                    green = 0;
                    blue = 0;
                  }
                  
                  labelCnt++;
                }
              }

              if(labelCnt >= labelWidth*labelHeight){
                triggerLabel = false;
              }
            }
            
            if(triggerLabel2 == true){
              
              if((labelYOffset2 <= y) && (y < labelYOffset2+labelHeight2)){
                UInt idRow = floor((y-labelYOffset2)/upscaleLabel);
              
                if((labelXOffset2 <= x) && (x < labelXOffset2+labelWidth2)){
                  UInt idCol = floor((x-labelXOffset2)/upscaleLabel);
//                  std::cout << "idCol, idRow: " << idCol << ", " << idRow << std::endl;
                  
                  // go over bitmap; 1 = black pixel
                  if(currentLabel2[idRow][idCol] == 1){
                    red = 0;
                    green = 0;
                    blue = 0;
                  }
                  
                  labelCnt2++;
                }
              }

              if(labelCnt2 >= labelWidth2*labelHeight2){
                triggerLabel2 = false;
              }
            }
            
                            
          }
        } else {
          /*
           * rotation state set
           * > coloring according to angle between rotation direction and x-axis (only x-y-plane considered!)
           * > color map over angle:
           *         0  - 120 -> red decreases linearly, blue increases linearly, green = 0
           *        120 - 240 -> red = 0, blue decreases linearly, green increases linearly
           *        240 - 360 -> red increases linearly, blue = 0, green decreases linearly
           * > switching state is encoded in brightness information
           * >      +1  -> full brightness, i.e., 255
           * >      -1  -> full brightness but every second pixel set to black
           */
          
          // angle between 0 and 360 degree needed > add 180
          // Note: atan2 takes y and x as argument 1 and 2 and computes the angle towards the x-axis
          //      > 
          //      
          Double angle = 0; 
          Double magnitude = 0;
          
          if(rotStateGiven){
            angle = std::atan2(rot_y,rot_x)/M_PI * 180 + 180;
            // put info of data into magnitude
            magnitude = scaling*255.0;
          } else {
            /*
             * color by normalized data_[idy][idx]/max(max(abs(data_))); 
             * to reuse algorithm, map -1.0 to 0 degree and 1.0 to 360 degree
             */
            angle = data_[idy][idx]/data_abs_max*180.0 + 180.0;
            // here info of data is already in angle > keep magnitude fix
            magnitude = 255.0;
          }
          
          GetColor_from_ColorMap(colormapId, angle, magnitude, &red, &green, &blue);          
//          
//          if(angle >= 0 && angle < 120){
//            red = lround( magnitude*(120-angle)/120.0 ) ;
//            blue = lround( magnitude*(angle/120.0 ) );
//            green = 0;
//          } else if(angle >= 120 && angle < 240){
//            red = 0;
//            blue = lround( magnitude*(240-angle)/120.0 );
//            green = lround( magnitude*(angle-120)/120.0 );
//          } else {
//            red = lround( magnitude*(angle-240)/120.0 );
//            blue = 0;
//            green = lround( magnitude*(360-angle)/120.0 );
//          }
//          
          if(data_[idy][idx] < 0){
            // every second pixel black
            if((x+y)%2 == 0){
              red = 0;
              blue = 0;
              green = 0;
            } 
          }
        }
        
        /*
         * print grid by drawing black lines of lineThickness pixel height
         */
        if(printGrid){
          UInt lineThickness = 2; // should be a multiple of 2!
          UInt extendOverTriangle = 15;
          Double numLines = 10.0;
          
          /*
           * horizontal lines
           */
          // at -1.0
          // only half linethickness - this is why lineThickness should be a multiple of 2
          if( (y >= 0) && (y <= ceil(lineThickness/2.0)) ){
            if( x <= y+extendOverTriangle ){
              red = 0;
              green = 0;
              blue = 0;
            }
          }
          
          for(UInt horizontalIdx = 1; horizontalIdx < numLines; horizontalIdx++){
          
            // at -0.8, -0.6, ... , 0.6, 0.8
            if( (y >= ceil(horizontalIdx*height/numLines)-lineThickness/2.0) && (y <= ceil(horizontalIdx*height/numLines)+lineThickness/2.0) ){
              if( x <= y+extendOverTriangle ){
                red = 0;
                green = 0;
                blue = 0;
              }
            }
          }
          
          // at 1.0
          if( (y >= height-ceil(lineThickness/2.0)) && (y <= height) ){
            if( x <= y+extendOverTriangle ){
              red = 0;
              green = 0;
              blue = 0;
            }
          }
        
          /*
           * vertical lines
           */
          // at -1.0
          // only half linethickness - this is why lineThickness should be a multiple of 2
          if( (x >= 0) && (x <= ceil(lineThickness/2.0)) ){
            if( y >= x-extendOverTriangle ){
              red = 0;
              green = 0;
              blue = 0;
            }
          }
          
          for(UInt verticalIdx = 1; verticalIdx < numLines; verticalIdx++){
          
            // at -0.8, -0.6, ... , 0.6, 0.8
            if( (x >= ceil(verticalIdx*width/numLines)-lineThickness/2.0) && (x <= ceil(verticalIdx*width/numLines)+lineThickness/2.0) ){
              if( y >= x-extendOverTriangle ){
                red = 0;
                green = 0;
                blue = 0;
              }
            }
          }
          
          // at 1.0
          if( (x >= width-ceil(lineThickness/2.0)) && (x <= width) ){
            if( y >= x-extendOverTriangle ){
              red = 0;
              green = 0;
              blue = 0;
            }
          }
          
          /*
           * add labeled square block indicating the blocksize/cellsize in the grid
           * > only draw contour
           */
          if( (y >= blockYOffset)&&(y < blockYOffset + blockHeight) ){
            if( ( (x >= blockXOffset-ceil(lineThickness/2.0))&&(x <= blockXOffset+ceil(lineThickness/2.0)) )  ||
                  ( (x >= blockXOffset+blockWidth-1-ceil(lineThickness/2.0))&&(x <= blockXOffset+blockWidth-1+ceil(lineThickness/2.0)) ) ) { 
              red = 0;
              green = 0;
              blue = 0;
            }
          }
          if( ( (y >= blockYOffset-ceil(lineThickness/2.0))&&(y <= blockYOffset+ceil(lineThickness/2.0)) )  ||
                ( (y >= blockYOffset+blockHeight-1-ceil(lineThickness/2.0))&&(y <= blockYOffset+blockHeight-1+ceil(lineThickness/2.0)) ) ) { 
            if( (x >= blockXOffset)&&(x < blockXOffset + blockWidth) ){
              red = 0;
              green = 0;
              blue = 0;
            }
          }
          
          // label on side
          if( ((y >= blockYOffset)&&(y <= blockYOffset+1)) && (triggerLabel3 == false) ) {
            triggerLabel3 = true;
            labelCnt3 = 0;
            currentLabel3 = zeroDotTwoSide;
            labelYOffset3 = blockYOffset + ceil(lineThickness/2.0);
            labelXOffset3 = blockXOffset + blockWidth + lineThickness;
          }
          
          // label on top
          if( ((y >= blockYOffset+blockHeight)&&(y <= blockYOffset+blockHeight+1)) && (triggerLabel3 == false) ) {
            triggerLabel3 = true;
            labelCnt3 = 0;
            currentLabel3 = zeroDotTwoTop;
            labelYOffset3 = blockYOffset + blockHeight + lineThickness;
            labelXOffset3 = blockXOffset + ceil(lineThickness/2.0);
          }
//              bool triggerLabel3 = false;
//    UInt labelCnt3 = 0;
//    static UInt labelBaseWidth3 = 12;
//    static UInt labelBaseHeight3 = 12;
//    
//    UInt labelWidth3 = labelBaseWidth3*upscaleLabel;
//    UInt labelHeight3 = labelBaseHeight3*upscaleLabel;
//    UInt blockWidth = ceil(width/20);
//    UInt blockHeight = blockWidth;        
//    UInt blockXOffset = ceil(0.4*colorBarXOffset);
//    UInt blockYOffset = colorBarYOffset;
//    UInt labelXOffset3 = blockXOffset + blockWidth;
//    UInt labelYOffset3 = blockYOffset + blockHeight;
          
            if(triggerLabel3 == true){
              
              if((labelYOffset3 <= y) && (y < labelYOffset3+labelHeight3)){
                UInt idRow = floor((y-labelYOffset3)/upscaleLabel);
              
                if((labelXOffset3 <= x) && (x < labelXOffset3+labelWidth3)){
                  UInt idCol = floor((x-labelXOffset3)/upscaleLabel);
//                  std::cout << "idCol, idRow: " << idCol << ", " << idRow << std::endl;
                  
                  // go over bitmap; 1 = black pixel
                  if(currentLabel3[idRow][idCol] == 1){
                    red = 0;
                    green = 0;
                    blue = 0;
                  }
                  
                  labelCnt3++;
                }
              }

              if(labelCnt3 >= labelWidth3*labelHeight3){
                triggerLabel3 = false;
              }
            }      
        }
          
//        if(idx == idy){
//          if(y < x){
//            /*
//             * make pixel white in lower half of main diagonal to form triangles
//             */
//            red = 255;
//            blue = 255;
//            green = 255;
//          }
//        }

        unsigned char pixel[3];
        pixel[0] = blue;
        pixel[1] = green;
        pixel[2] = red;

        outfile.write( (char*)pixel, 3 );
      }
      outfile.write( (char*)pad, padsize );
    }

    outfile.close();
    
    delete[] colorbarTitle;
    delete[] unsetRotState;
    delete[] bmpNeg180;
    delete[] bmpNeg90;
    delete[] bmp0;
    delete[] bmpPos90;
    delete[] bmpPos180;
    
    delete[] switchingTitle;
    delete[] plus1;
    delete[] minus1;
    
    delete[] zeroDotTwoTop;
    delete[] zeroDotTwoSide;
  }

  template<>
  void Matrix<Integer>::PerformRotation( const Matrix<Double>& R,  Matrix<Integer>& retMat ) const {
    EXCEPTION("Rotation only defined for double- and complex valued matrixes");
  }

  template<>
  void Matrix<UInt>::PerformRotation( const Matrix<Double>& R,  Matrix<UInt>& retMat ) const {
    EXCEPTION("Rotation only defined for double- and complex valued matrixes");
  }
  
  template<>
  void Matrix<bool>::PerformRotation( const Matrix<Double>& R,  Matrix<bool>& retMat ) const {
    EXCEPTION("Rotation only defined for double- and complex valued matrixes");
  }

  
  template<class TYPE>
  void Matrix<TYPE>::PerformRotation( const Matrix<Double>& R,  Matrix<TYPE>& retMat ) const {
    
    // Note; Currently the rotation only works for 3x3, 3x6, 6x3 and 6x6 matrices.
    // However, we should generalize the rotation also for 2x4 and 4x4 matrices for the
    // 2D and axi case for mechanics.

    //get dimension of matrix
    UInt rowSize = size_row_;
    UInt colSize = size_col_;

    Matrix<TYPE> helpMat;

    if ( rowSize == 3 && colSize == 3 && R.GetNumCols() == 3 && R.GetNumRows() == 3) {
      // get memory for transposed rotation matrix
      Matrix<Double> RT;
      RT.Resize(3,3);
      R.Transpose(RT);
      // tensor is a 3x3 matrix: sol = R * matrixOrig * RT
      helpMat   = (*this) * RT;
      retMat = R * helpMat;
    } else if ( rowSize == 2 && colSize == 2 && R.GetNumCols() == 2 && R.GetNumRows() == 2) {
      // get memory for transposed rotation matrix
      Matrix<Double> RT;
      RT.Resize(2,2);
      R.Transpose(RT);
      // tensor is a 3x3 matrix: sol = R * matrixOrig * RT
      helpMat   = (*this) * RT;
      retMat = R * helpMat;
      
    } else if(((rowSize == 2 && colSize == 3) ||
               (rowSize == 3 && colSize == 2) ||
               (rowSize == 3 && colSize == 3)) &&
             (R.GetNumCols() == 2 && R.GetNumRows() == 2)) {
      // case of 2x3 coupling tensors as well as 3x3 tensors with a 2x2 rotation matrix
      // 2D tensor rotation
      Matrix<Double> Q;
      Q.Resize(3,3);
      Q[0][0] = R[0][0]*R[0][0];
      Q[0][1] = R[0][1]*R[0][1];
      Q[0][2] = 2.0*R[0][0]*R[0][1];
      Q[1][0] = R[1][0]*R[1][0];
      Q[1][1] = R[1][1]*R[1][1];
      Q[1][2] = 2.0*R[1][0]*R[1][1];
      Q[2][0] = R[0][0]*R[1][0];
      Q[2][1] = R[0][1]*R[1][1];
      Q[2][2] = R[0][0]*R[1][1] + R[0][1]*R[1][0];

      Matrix<Double> QT;
      QT.Resize(3,3);
      Q.Transpose(QT);

      // shouldn't this work for 2x3 tensors?
      if ( rowSize == 2 && colSize == 3 ) {
        helpMat   = (*this) * QT;
        retMat = R * helpMat;
      }
      else if (rowSize == 3 && colSize == 3 ) {
        helpMat   = (*this) * QT;
        retMat = Q * helpMat;
      }
      else if ( rowSize == 2 && colSize == 3 ) {
        Matrix<Double> RT;
        RT.Resize(2,2);
        R.Transpose(RT);
        helpMat   = (*this) * RT;
        retMat = Q * helpMat;
      }

      //helpMat   = (*this) * QT;
      //retMat = Q * helpMat;
//      std::cout<<"Q = "<<Q.ToString()<<std::endl;

    } else if( (rowSize == 3 && colSize == 6) ||
             (rowSize == 6 && rowSize == 6 ) ) {
      // we also need Q;
      Matrix<Double> Q;

      // Composed Rotation Matrix
      // Ref.: M.Richter, "Entwicklung mechanischer Modelle zur analytischen
      // Beschreibung der Materialeigenschaften von textilbewehrtem Feinbeton",
      // Diss., Dresden, 2005, p. 27

      Q.Resize(6,6);  

      Q[0][0] = R[0][0]*R[0][0];
      Q[0][1] = R[0][1]*R[0][1];
      Q[0][2] = R[0][2]*R[0][2];
      Q[0][3] = 2.0*R[0][1]*R[0][2];
      Q[0][4] = 2.0*R[0][0]*R[0][2];
      Q[0][5] = 2.0*R[0][0]*R[0][1];

      Q[1][0] = R[1][0]*R[1][0];
      Q[1][1] = R[1][1]*R[1][1];
      Q[1][2] = R[1][2]*R[1][2];
      Q[1][3] = 2.0*R[1][1]*R[1][2];
      Q[1][4] = 2.0*R[1][0]*R[1][2];
      Q[1][5] = 2.0*R[1][0]*R[1][1];

      Q[2][0] = R[2][0]*R[2][0];
      Q[2][1] = R[2][1]*R[2][1];
      Q[2][2] = R[2][2]*R[2][2];
      Q[2][3] = 2.0*R[2][1]*R[2][2];
      Q[2][4] = 2.0*R[2][0]*R[2][2];
      Q[2][5] = 2.0*R[2][0]*R[2][1];

      Q[3][0] = R[1][0]*R[2][0];
      Q[3][1] = R[1][1]*R[2][1];
      Q[3][2] = R[1][2]*R[2][2];
      Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
      Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
      Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];

      Q[4][0] = R[0][0]*R[2][0];
      Q[4][1] = R[0][1]*R[2][1];
      Q[4][2] = R[0][2]*R[2][2];
      Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
      Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
      Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];

      Q[5][0] = R[0][0]*R[1][0];
      Q[5][1] = R[0][1]*R[1][1];
      Q[5][2] = R[0][2]*R[1][2];
      Q[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
      Q[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0];
      Q[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0];


      //  std::cout << "R:\n" << R << std::endl;
      //  std::cout << "Q:\n" << Q << std::endl;
      //  std::cout << "Tensor orig:\n" << matTensor << std::endl;

      Matrix<Double> QT;
      QT.Resize(6,6);
      Q.Transpose(QT);

      if ( rowSize == 3 && colSize == 6 ) {
        helpMat   = (*this) * QT;
        retMat = R * helpMat;
      }
      else if (rowSize == 6 && colSize == 6 ) {
        helpMat   = (*this) * QT;
        retMat = Q * helpMat;
      }
    }
    else if (rowSize == 1 && colSize == 1) {
      retMat = *this; // No rotation needed
    }
    else {
      EXCEPTION("Rotation of " << rowSize << " x " << colSize << " matrix not implemented yet.");
    }

  }
  
  template<class TYPE>
  unsigned int Matrix<TYPE>::ParseLineHelper(const std::string& input, StdVector<TYPE>& out)
  {
    // the rows are in the form "1 2 3] : " or "1.2 2.5 3.5] : " or "(1,1) (0,-2)] : "

    // out will be resized in all but the first calls. So Push_back is not that expensive!
    char_separator<char> sep(" ]:"); // spaces, closing and colon
    tokenizer<char_separator<char> > tokens(input, sep);

    for(tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it)
    {
      try
      {
        out.Push_back(boost::lexical_cast<TYPE>(*it));
      }
      catch(boost::bad_lexical_cast &)
      {
        EXCEPTION("Cannot cast to value while parsing matrix: '" << *it << "'");
      }
    }

    if(out.GetSize() == 0) EXCEPTION("cannot interpret matrix row as data: '" << input << "'");

    return out.GetSize();
  }


  template<class TYPE>
  void Matrix<TYPE>::Parse(const std::string& input)
  {
    // we don't know the number or rows and cols in advance.
    StdVector<std::string> rows;
    char_separator<char> sep("["); // ignore the row number and colon
    tokenizer<char_separator<char> > tokens(input, sep);

    for(tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it)
      rows.Push_back(*it);

    // the rows are in the form "1 2 3] : " or "1.2 2.5 3.5] : " or "(1,1) (0,-2)] : "
    if(rows.GetSize() == 0) EXCEPTION("cannot interpret data as matrix with data: '" << input << "'");
    StdVector<TYPE> row;
    // check first row
    unsigned int n_cols = ParseLineHelper(rows[0], row); // row is now set to the proper size

    // we have the information, store the data!
    this->Resize(rows.GetSize(), n_cols);

    for(unsigned int r = 0; r < rows.GetSize(); r++)
    {
      row.Clear();
      row.Reserve(n_cols); // I know no Clear() Reserve() combination, but it's not too expensive

      ParseLineHelper(rows[r], row); // row is now set to the proper size
      if(row.GetSize() != n_cols)
        EXCEPTION("matrix has inconsistent number or columns within row " << (r+1) << ": '" << input << "'");

      for(unsigned int c = 0; c < row.GetSize(); c++)
        data_[r][c] = row[c];
    }
  }


  template<class TYPE>
  void Matrix<TYPE>::Resize(const UInt nRows, const UInt nCols )
  {
    if(nRows != size_row_ || nCols != size_col_)
    {
      // delete old data, if existend
      if (data_ != NULL)
      {
        delete [] data_[0];
        delete [] data_;
      }

      // set the size to requested values
      size_row_ = nRows; 
      size_col_ = nCols;
      
      if (size_row_ == 0 || size_col_ == 0)
      {
        data_ = nullptr;
      }
      else {
        data_    = new TYPE*[size_row_];
        data_[0] = new TYPE [size_row_*size_col_];

        for (UInt k = 1; k < size_row_; ++k) 
          data_[k] = data_[k-1] + size_col_;

      }
    }
  }


  template<class TYPE>
  void Matrix<TYPE>::Resize(const UInt col )
  {
    Resize(col,col);  
  }

  template<class TYPE>
  void Matrix<TYPE>::Resize(const Matrix<TYPE>& other)
  {
    Resize(other.size_row_, other.size_col_);  
  }

#ifndef USE_EXPRESSION_TEMPLATES

  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator=(const Matrix<TYPE> &x)
  {
    if (this == &x)
      return *this;

    // set the size in any case to the size of the assigned matrix
    Resize(x.size_row_, x.size_col_);
  
    // allows to copy an empty matrix. !
    if (size_row_ > 0 && size_col_ > 0) {
      // copy the entries
      std::memcpy(data_[0], x.data_[0], size_row_ * size_col_ * sizeof(TYPE));
    }

    return *this;
  }

  template<class TYPE>
  Matrix<TYPE> Matrix<TYPE>::operator+ () const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif

    return *this;
  }

 

  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator+=(const Matrix<TYPE> &x)
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif
  
#ifdef CHECK_INDEX  
    if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
      EXCEPTION("incompatible dimension for +"); 
#endif
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] += x.data_[0][k];
  
    return *this;
  }

  template<class TYPE>
  Matrix<TYPE> Matrix<TYPE>::operator-() const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 )
      EXCEPTION("undefined Matrix");
#endif

    Matrix<TYPE> z(size_row_,size_col_);
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        z[0][k] = -data_[0][k];
    
    return z;
  }

 

  template<class TYPE>
  Matrix<TYPE> & Matrix<TYPE>::operator-=(const Matrix<TYPE> &x)
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif
  
#ifdef CHECK_INDEX  
    if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
      EXCEPTION("incompatible dimension for +"); 
#endif
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] -= x.data_[0][k];
    
    return *this;
  }


  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator/= (const TYPE x)
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k] /= x;

    return *this;
  }


#endif // USE_EXPRESSION_TEMPLATES
  
  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator*= (const TYPE x)
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif
    UInt size = size_row_ * size_col_;
    for(UInt k = 0, s = size; k < s; ++k)
      data_[0][k] *= x;

    return *this;
  }

 template<class TYPE>
  Matrix<TYPE> & Matrix<TYPE>::operator*=(const Matrix<TYPE> &x)
  {   
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX  
    if (size_col_ != x.size_row_)
      EXCEPTION("incompatible dimension");
#endif
 
    TYPE    a;
    Matrix  z (size_row_, x.size_col_);
  
    UInt i,j; 
    for (i = 0; i < size_row_; i++)
      for (j = 0; j < x.size_col_; j++)
        {       
          a = data_ [i] [0] * x.data_ [0] [j];
          for (UInt k = 1; k < size_col_; k++)
            a += data_ [i] [k] * x.data_ [k] [j];
          z.data_ [i] [j] = a;
        }
  
    *this = z;
    return *this;
  }

  template<class TYPE>
  bool Matrix<TYPE>::operator== (const Matrix<TYPE> &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if(data_ [0][k] != x.data_[0][k]) return false;
  
    return true;
  }

  template<class TYPE>
  bool Matrix<TYPE>::operator!= (const Matrix<TYPE> &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif 
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if (data_ [0][k] != x.data_[0][k]) return true;
  
    return false;
  }
  
  /** Adds the multiple of another matrix */
  template<class TYPE>
  void Matrix<TYPE>::Add(const TYPE factor, const Matrix<TYPE>& mat)
  {
    const Matrix<TYPE>& other_mat = dynamic_cast<const Matrix<TYPE> & >(mat);
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || other_mat.size_row_ == 0 || other_mat.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX
    if (size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_)
      EXCEPTION("matrices do not match");
#endif
    
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] += factor * other_mat.data_[0][k];
  }

  /** Adds the multiple of the transpose of another matrix */
  template<class TYPE>
  void Matrix<TYPE>::AddT(const TYPE factor, const Matrix<TYPE>& mat)
  {
    const Matrix<TYPE>& other_mat = dynamic_cast<const Matrix<TYPE> & >(mat);
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || other_mat.size_row_ == 0 || other_mat.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX
    if (size_row_ != other_mat.size_col_ || size_col_ != other_mat.size_row_)
      EXCEPTION("matrices do not match");
#endif

    for(UInt k = 0, s = size_row_ ; k < s; ++k)
      for (UInt kk =0, ss = size_col_; kk < ss; ++kk)
        data_[k][kk] += factor * other_mat.data_[kk][k];
  }

  /** Assigns a multiple of another matrix */
  template<class TYPE>
  void Matrix<TYPE>::Assign(const Matrix<TYPE>& other_mat, TYPE factor, bool size_tolerant)
  {
    if(size_row_ == other_mat.size_row_ && size_col_ == other_mat.size_col_)
      for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] = factor * other_mat.data_[0][k];
    else
    {
      if(!size_tolerant)
        EXCEPTION("matrices do not match");
      // assure the non assigned area is zero
      if(size_row_ > other_mat.size_row_ || size_col_ > other_mat.size_col_)
        this->Init();
      for(unsigned int r = 0, nr = std::min(size_row_, other_mat.size_row_); r < nr; r++)
        for(unsigned int c = 0, nc = std::min(size_col_, other_mat.size_col_); c < nc; c++)
          data_[r][c] = factor * other_mat.data_[r][c];
    }
  }
  
  template<class TYPE>
  void Matrix<TYPE>::Assign(const Vector<TYPE>& vec, unsigned int rows, unsigned int cols, bool row_major)
  {
    assert(vec.GetSize() == rows * cols);
    assert(vec.GetSize() > 0);

    Resize(rows, cols);

    for(unsigned int r = 0; r < rows; r++)
      for(unsigned int c = 0; c < cols; c++)
        data_[r][c] = row_major ? vec[r*cols + c] : vec[c*rows+r];
  }
  
  // perform matrix-matrix multiplication using BLAS (general case)
  template<class TYPE>
  void Matrix<TYPE>::Mult_Blas( const Matrix<TYPE>& mMat1,
                                Matrix<TYPE>& rMat1,
                                bool trans_a, bool trans_b,
                                TYPE alpha, TYPE beta, bool conjugate ) const {

#ifdef CHECK_INDEX
    if((trans_a == true) && (trans_b == true)){
      if (size_row_ != mMat1.GetNumCols())
        EXCEPTION("incompatible dimension");
    } else if((trans_a == false) && (trans_b == true)){
      if (size_col_ != mMat1.GetNumCols())
        EXCEPTION("incompatible dimension");
    } else if((trans_a == true) && (trans_b == false)){
      if (size_row_ != mMat1.GetNumRows())
        EXCEPTION("incompatible dimension");
    } else {
      if (size_col_ != mMat1.GetNumRows())
        EXCEPTION("incompatible dimension");
    }
#endif

#ifdef CHECK_INITIALIZED
    if(size_row_ == 0 || size_col_ == 0)
      EXCEPTION("undefined Matrix");
    if(mMat1.GetNumRows() == 0 || mMat1.GetNumCols() ==0)
      EXCEPTION("undefined Matrix");
    if(rMat1.GetNumRows() == 0||rMat1.GetNumCols()==0)
      EXCEPTION("undefined Matrix");
#endif

    // because of the column-wise access of fortran the routine would calc
    // C^T = A^T * B^T if you give it C, A and B
    // but because A^T * B^T = (B * A)^T you can just calculate:
    // C = B * A
    // --> swap A and B

    // actually many transposed in complex are meant to be conjugate complex, e.g. B'(DB). Unfortunately
    // this is often not considered in CFS, so take care.
    char transp = conjugate && boost::is_complex<TYPE>::value ? 'c' : 't';

    char transa = trans_a ? transp : 'n';
    char transb = trans_b ? transp : 'n';

    // m would normally be the number of rows of C but Fortran accesses it as C^T so m is the number of columns
    //  in the same way n is now the number of rows
    //  k is the number of rows of op(B) and in our case op(B) = A or A^T so k would be rows of A in the first case and cols of A in the second one
    // but here you also have to remember the column wise access so cols and rows are swapped like for c

    int n = rMat1.GetNumRows();
    int m = rMat1.GetNumCols();
    int k = trans_a ? size_row_ : size_col_;

    TYPE* A = data_[0];
    TYPE* B = mMat1.data_[0];
    TYPE* C = rMat1.data_[0];

    // here you use the same properties as in the documentation of dgemm with the difference,
    // that our lda is LDB and ldb is LDA because we swapped A and B

    int lda = trans_a ? n : k;
    int ldb = trans_b ? k : m;
    int ldc = m;
    CallGEMM(&transb,&transa,&m,&n,&k,&alpha,B,&ldb,A,&lda,&beta,C,&ldc);
   }
  

  template<class TYPE >
  void Matrix<TYPE>::CallGEMM(char* transb, char* transa, int* m, int* n, int* k, TYPE* alpha, TYPE* B, int* ldb, TYPE* A, int* lda, TYPE* beta, TYPE* C, int* ldc) const {
    EXCEPTION("wrong type");
  }

  template< >
  void Matrix<Double>::CallGEMM(char* transb, char* transa, int* m, int* n, int* k, Double* alpha, Double* B, int* ldb, Double* A, int* lda, Double* beta, Double* C, int* ldc) const {
    dgemm(transb,transa,m,n,k,alpha,B,ldb,A,lda,beta,C,ldc);
  }

  template< >
  void Matrix<Complex>::CallGEMM(char* transb, char* transa, int* m, int* n, int* k, Complex* alpha, Complex* B, int* ldb, Complex* A, int* lda, Complex* beta, Complex* C, int* ldc) const {
    zgemm(transb,transa,m,n,k,alpha,B,ldb,A,lda,beta,C,ldc);
  }

  // Perform a matrix-vector multiplication rvec = this*mvec
  template<class TYPE>
  void Matrix<TYPE>::Mult(const SingleVector & mvec, SingleVector & rvec) const
  {
    Vector<TYPE> const & mvec1 = dynamic_cast<const Vector<TYPE>& >(mvec);
    Vector<TYPE> & rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
#if defined CHECK_INITIALIZED || defined CHECK_INDEX
    UInt size_mvec = mvec1.GetSize();
    UInt size_rvec = rvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mvec == 0) 
      EXCEPTION("undefined Vector");
    if (size_rvec == 0) 
      EXCEPTION("undefined Vector");
#endif

#ifdef CHECK_INDEX
    if (size_col_ != size_mvec) 
      EXCEPTION("incompatible dimension");
    if (size_row_ != size_rvec) 
      EXCEPTION("incompatible dimension");
#endif

#endif
    
    UInt k,kk;
    rvec1.Init();
    for ( k = 0; k < size_row_; k++)
      for ( kk = 0; kk < size_col_; kk++)
        rvec1[k] += data_[k][kk]*mvec1[kk];
  }


  template<class TYPE>
  void Matrix<TYPE>::Mult(TYPE scale)
  {
    for(unsigned int i = 0; i < size_row_ * size_col_; i++)
      data_[0][i] *= scale;
  }


  template<class TYPE>
  TYPE Matrix<TYPE>::ScalarProduct(const Matrix<TYPE>& other_mat) const
  {
#ifdef CHECK_INITIALIZED
    if(size_row_ == 0 || size_col_ == 0)
      EXCEPTION("undefined Matrix");
    if(size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_)
      EXCEPTION("incompatible dimension");
#endif

    TYPE result(0);

    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      result += OpType<TYPE>::dotProduct(data_[0][k], other_mat.data_[0][k]);

    return result;
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::FrobeniusProduct(const Matrix<TYPE>& other_mat) const
  {
#ifdef CHECK_INITIALIZED
    if(size_row_ == 0 || size_col_ == 0)
      EXCEPTION("undefined Matrix");
    if(size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_)
      EXCEPTION("incompatible dimension");
#endif

    TYPE result(0);

    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      result += OpType<TYPE>::dotProduct(data_[0][k], other_mat.data_[0][k]);

    return result;
  }

  template<>
  Matrix<double> Matrix<double>::EntryMult(const Matrix<double>& other_mat) const
  {
#ifdef CHECK_INITIALIZED
    if(size_row_ == 0 || size_col_ == 0)
      EXCEPTION("undefined Matrix");
    if(size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_)
      EXCEPTION("incompatible dimension");
#endif

    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k] *= other_mat.data_[0][k];

    return *this;
  }



  // Perform a matrix-vector multiplication rvec = this*mvec via the inner product
  template<class TYPE>
  void Matrix<TYPE>::MultInner(const SingleVector & mvec, SingleVector & rvec) const
  {
    Vector<TYPE> const & mvec1 = dynamic_cast<const Vector<TYPE>& >(mvec);
    Vector<TYPE> & rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
#if defined CHECK_INITIALIZED || defined CHECK_INDEX
    UInt size_mvec = mvec1.GetSize();
    UInt size_rvec = rvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mvec == 0) 
      EXCEPTION("undefined Vector");
    if (size_rvec == 0) 
      EXCEPTION("undefined Vector");
#endif

#ifdef CHECK_INDEX
    if (size_col_ != size_mvec) 
      EXCEPTION("incompatible dimension");
    if (size_row_ != size_rvec) 
      EXCEPTION("incompatible dimension");
#endif

#endif

    UInt k,kk;
    rvec1.Init();
    for ( k = 0; k < size_row_; k++)
      for ( kk = 0; kk < size_col_; kk++)
        rvec1[k] += OpType<TYPE>::dotProduct(data_[k][kk],mvec1[kk]);
  }
  
  // Perform a matrix-vector multiplication rvec = transpose(this)*mvec
  template<class TYPE>
  void Matrix<TYPE>::MultT(const SingleVector &mvec, SingleVector &rvec) const
  {
    Vector<TYPE> const &mvec1 = dynamic_cast<const Vector<TYPE>& >(mvec);
    Vector<TYPE> &rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
#if defined CHECK_INITIALIZED || defined CHECK_INDEX
    UInt size_mvec = mvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mvec == 0) 
      EXCEPTION("undefined input Vector");
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_mvec) 
      EXCEPTION("incompatible dimension");
#endif

#endif
    
    // overwrite output vector with 0.0s and set correct length
    rvec1.Resize(size_col_, 0);
    for ( UInt k = 0; k < size_row_; ++k)
      for ( UInt kk = 0; kk < size_col_; ++kk)
        rvec1[kk] += data_[k][kk]*mvec1[k];
  }
  
  // **************
  //   operator<<
  // **************
  template<class S>
  std::ostream &operator<< ( std::ostream &out, const Matrix<S> &mat ) {


    // Set output format
    out.setf( std::ios::scientific );
    const UInt oldPrec = out.precision();
    const UInt newPrec = 7;
    out.precision( newPrec );

    const unsigned int mrows(mat.GetNumRows());
    const unsigned int mcols(mat.GetNumCols());
    
    for ( UInt i = 0; i < mrows; ++i ) {
      for ( UInt j = 0; j < mcols; ++j ) {
        out << std::setw( newPrec + 7 ) << mat[i][j] << " ";
      }
      out << std::endl;
    }

    // Restore old settings
    out.precision( oldPrec );

    return out;
  }

  template std::ostream &operator<<<UInt>    ( std::ostream &,
                                               const Matrix<UInt>    & );
  template std::ostream &operator<<<Integer> ( std::ostream &,
                                               const Matrix<Integer> & );
  template std::ostream &operator<<<Double>  ( std::ostream &,
                                               const Matrix<Double>  & );
  template std::ostream &operator<<<Complex> ( std::ostream &,
                                               const Matrix<Complex> & );


  template<class TYPE>
  void Matrix<TYPE>::Transpose (Matrix<TYPE> &transposedMat) const
  {
    transposedMat.Resize(size_col_, size_row_);
 
    UInt i,j;
 
    for (i = 0; i < size_col_; i++)
      for (j = 0; j < size_row_; j++)
        transposedMat.data_[i][j] = data_[j][i];
  }

  

  template<class TYPE>
  bool Matrix<TYPE>::IsHermitian(double eps) const
  {
    for(UInt r = 0; r < size_row_; r++)
      for(UInt c = r+1; c < size_col_; c++)
        if(!close(data_[r][c], data_[c][r], eps))
          return false;

    return true;
  }


  template< >
  bool Matrix<Complex>::IsHermitian(double eps) const
  {
    for(UInt r = 0; r < size_row_; r++)
    {
      if(!IsNoise(data_[r][r].imag()))
        return false;

      for(UInt c = r+1; c < size_col_; c++)
      {
        Complex u = data_[r][c]; // upper
        Complex l = data_[c][r]; // lower

        if(!close(u.real(), l.real()))
          return false;

        if(!close(std::abs(u.imag()), std::abs(l.imag())))
          return false;
      }
    }
    return true;
  }


  template<class TYPE>
  bool Matrix<TYPE>::DirectSolve(SingleVector& x1, const SingleVector& b1, bool throw_exception)
  {
    assert(!ContainsInf() && !ContainsNaN());

    Vector<TYPE> & x = dynamic_cast<Vector<TYPE>& >(x1);
    const Vector<TYPE> & b = dynamic_cast<const Vector<TYPE>& >(b1);

    x.Resize(size_row_);

    Integer nmat = size_row_-1;
    Integer i, j, k, k1;

    //  the Gauss elimination 
    for(k = 0; k < nmat; ++k)
    {
      k1 = k + 1;
      for(i = k1; i <= nmat; ++i)
      {
        if (data_[k][k] != 0.0)
        {
          data_[i][k] /= data_[k][k];
          for (j=k1; j<=nmat; ++j)
            data_[i][j] -= data_[i][k] * data_[k][j];
        }
        else
        {
          if(throw_exception)
            throw Exception("Cannot invert zero element at step " + std::to_string(k));
          else
            return false;
        }
      }
    }

    // TODO: shall work w/o temporary vector. See e.g. petsc toptop implementation from Niels Aage
    // solve Ly = b by forward substitution
    Vector<TYPE> y(b.GetSize());

    for (i=0; i<=nmat; ++i)
    {
      y[i] = b[i];
      for(j = 0; j < i; ++j)
        y[i] -= data_[i][j] * y[j];
    }

    // solve Ux = y backward substitution
    
    for(i = nmat; i >= 0; --i)
    {
      x[i] = y[i];
      for(j = nmat; j > i; --j)
        x[i] -= data_[i][j] * x[j];
      x[i] /= data_[i][i];
    }

    return true;
  }


  template<>
  void Matrix<Complex>::solveWithLapack(Matrix<Complex> & b1,
                                        lapackSysMatType & LAPACK_MATRIX_TYPE)
  {

    if(size_row_!=size_col_)
      EXCEPTION("Matrix is not quadratic, no solver!" );

    char lp_matType;
    lp_matType='L';
    
    Integer lp_nrRHS, lp_info, lp_lwork,lp_dim;
    Integer lp_lda, lp_ldb;
    const unsigned int bcols(b1.GetNumCols());
    lp_nrRHS=bcols;
    lp_dim=size_row_;
    lp_lda=size_row_;
    lp_ldb=size_row_;

    Integer *lp_interchanges;
    
    std::complex<double>* lp_rhsVecf77;
    std::complex<double>* lp_sysVecf77;
    std::complex<double>* lp_workf77;
    
    Vector<Complex> lp_sysVec, lp_work;
    lp_sysVec.Resize(size_row_*size_row_);
    
      // copy values from system and RHS - Matrix into vector
    for (UInt i=0;i<size_row_;++i)
      for (UInt j=0;j<size_row_;++j){
        lp_sysVec[i+j*size_row_]=data_[i][j];
      }
    
    Vector<Complex> lp_rhsVec;
    const unsigned int brows(b1.GetNumRows());
    
    lp_rhsVec.Resize(brows*bcols);
 
    for (UInt i=0;i<brows;++i)
      for (UInt j=0;j<bcols;++j){
        lp_rhsVec[i+j*brows]=b1[i][j];
      }
    
    lp_rhsVecf77 = new std::complex<double>[size_row_*lp_nrRHS];
    lp_interchanges = new int[size_row_*size_row_];
    lp_workf77 = new std::complex<double>[size_row_*size_row_];
    lp_sysVecf77 = new std::complex<double>[size_row_*size_row_];
   

    // Convert openCFS Vector<Complex> to Vector<std::complex<double>>
    for ( UInt count = 0; count < size_row_*bcols; ++count ) {
      lp_rhsVecf77[count] = lp_rhsVec[count];
      }
    
    for (UInt count = 0; count < size_row_*size_row_; ++count ) {
      lp_sysVecf77[count] = lp_sysVec[count];
    }
    
    for (UInt count=0, ss = lp_work.GetSize(); count < ss; ++count ) {
      lp_workf77[count] = lp_work[count];
    }
    
    switch (LAPACK_MATRIX_TYPE){
      
    case ZGESVMTX:
      // solves systems with general system matrix
      zgesv(&lp_dim , &lp_nrRHS, lp_sysVecf77, &lp_lda, 
            lp_interchanges, lp_rhsVecf77, &lp_ldb, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZGESV reports invalid input parameter" );
      }    
      break;
    case ZSYSVMTX:
      
      lp_lwork=192;
      // solves systems with symmetric system matrix
      zsysv(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77, 
            &lp_lda, lp_interchanges, lp_rhsVecf77, &lp_ldb,
            lp_workf77, &lp_lwork, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZSYSV reports invalid input parameter" );
      }
      break;
    case ZHESVMTX:
      lp_lwork=192;
      // solves systems with hermitian system matrix
      zhesv(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77,
            &lp_lda, lp_interchanges,lp_rhsVecf77, &lp_ldb, 
            lp_workf77, &lp_lwork, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZHESV reports invalid input parameter" );
      }
      break;

    default:
      std::cout<<LAPACK_MATRIX_TYPE<<std::endl;
      EXCEPTION( " the chosen routine is not yet implemented in solveWithLapack" );


    } // matches switch ()...
      
          
     //reconvert Fortran77 -> CFS ++ datatypes
    for ( UInt count = 0; count < size_row_*bcols; ++count ) 
      lp_rhsVec[count] = lp_rhsVecf77[count];
          
    for ( UInt count = 0, ss = size_row_*size_row_; count < ss; ++count ) 
      lp_sysVec[count] = lp_sysVecf77[count];
          
    for (UInt count=0, ss = lp_work.GetSize(); count < ss; ++count )
      lp_work[count] = lp_workf77[count];
      
    // Writes result into b1
    for (UInt i=0;i<size_row_;++i)
      for (UInt j=0;j<bcols;++j)
        b1[i][j]=lp_rhsVec[i+j*brows];
  
    delete[] lp_rhsVecf77;
    delete[] lp_interchanges;
    delete[] lp_sysVecf77;
    delete[] lp_workf77;
    
  }


  template <class T>
  void Matrix<T>::eigenvaluesWithLapack(Vector<Double> & lp_w, Matrix<double> * ev_vec)
  {
    // computes all eigenvalues of a complex hermitian matrix

    char lp_jobz='V';
    char lp_uplo='L';
    Integer lp_N=size_row_;                
    Integer lp_lda=size_row_;
      
    // array contains ev in ascending order
    //    Vector<Double> lp_w;
    lp_w.Resize(size_row_);
    lp_w.Init();
    if (ev_vec != NULL) {
      (*ev_vec).Resize(size_row_,size_row_);
    }
    Integer lp_lworkf77=99;
      
    // workspace array - complex 16 array
    Vector<Complex> lp_work;
    lp_work.Resize(lp_lworkf77);
    lp_work.Init();
      
    // workspace array - double precision
    Vector<Double> lp_rwork;
    lp_rwork.Resize(3*size_row_-2);
    lp_rwork.Init();
      
    Integer lp_infof77;

    std::complex<double>* lp_af77    = new std::complex<double>[size_row_*size_row_];
    std::complex<double>* lp_workf77 = new std::complex<double>[99];
    double* lp_wf77     = new double[size_row_];
    double* lp_rworkf77 = new double[3*size_row_-2];
      
    // Convert openCFS Vector<Complex> to Vector<std::complex<double>>
    for ( UInt count = 0; count < size_row_; count++ ) 
      for ( UInt countC = 0; countC <size_row_; countC++ ) {
        lp_af77[countC*size_row_+count] = data_[count][countC];
      }
    
    for ( UInt count = 0; count < size_row_; count++ ) {
      lp_wf77[count] = lp_w[count];
    }
    
    for (UInt count=0; count < lp_rwork.GetSize();count++){
      lp_rworkf77[count] = lp_rwork[count];
    }
    
    zheev( &lp_jobz, &lp_uplo, &lp_N, lp_af77, &lp_lda, lp_wf77, 
           lp_workf77, &lp_lworkf77, lp_rworkf77 ,&lp_infof77); 
    
    // reconvert f772C++
    for (UInt count=0; count < lp_work.GetSize();count++)
      lp_work[count] = lp_workf77[count];
    
    for ( UInt count = 0; count < size_row_; count++ ) 
      lp_w[count] = lp_wf77[count];
    
    UInt c=0;
    if (ev_vec != NULL) {
      for ( UInt count = 0; count < size_row_; count++ ) {
        for (UInt count2 = 0; count2 < size_row_;count2++) {
          (*ev_vec)[count][count2] = lp_af77[c].real();
          if (std::abs(lp_af77[c].imag()) > std::numeric_limits<float>::epsilon() ) {
            EXCEPTION("Eigenvector is non real! ")
          }
          c++;
        }
      }
    }

    delete[] lp_workf77;
    delete[] lp_rworkf77;
    delete[] lp_af77;
    delete[] lp_wf77;
    
  }


  template<class TYPE>
  void Matrix<TYPE>::DyadicMult(const SingleVector & v1, const SingleVector & v2)
  {
    Vector<TYPE> const & vec1 = dynamic_cast<const Vector<TYPE>& >(v1);
    Vector<TYPE> const & vec2 = dynamic_cast<const Vector<TYPE>& >(v2);
  
    const UInt row = vec1.GetSize();
    const UInt col = vec2.GetSize();
  
    Resize(row, col);

    for(UInt actRow=0; actRow<row; actRow++)
      for(UInt actCol=0; actCol < col; actCol++)
        data_[actRow][actCol] = vec1[actRow] * vec2[actCol];
  }

  template<class TYPE>
  void Matrix<TYPE>::PseudoInvert(Matrix <TYPE> & inv) const
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0)
      EXCEPTION( "Undefined Matrix!" );
#endif

    Matrix<TYPE> metric;
    Matrix<TYPE> metricInverse;
    Matrix<TYPE> tmp;

    //std::cout<<*this<<std::endl;
    metric.Resize(size_col_, size_col_);
    MultT(*this, metric);
    //std::cout<<metric<<std::endl;
    metric.Invert(metricInverse);
    //std::cout<<metricInverse<<std::endl;
    tmp.Resize(size_row_, size_col_);
    Mult(metricInverse,tmp);
    tmp.Transpose(inv);
    //std::cout<<inv<<std::endl;

  }

  template<class TYPE>
  void Matrix<TYPE>::Invert (Matrix <TYPE> & inv) const
  {   

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION( "Undefined Matrix!" );
#endif

    if(size_row_== size_col_) {
      TYPE det;
      TYPE invDet;
      switch (size_row_)
        {
        case 1:
          inv.Resize(1);
          inv[0][0] = 1.0 / data_[0][0];
          break;
        case 2:
          inv.Resize(2,2);
          inv[0][0] = data_[1][1];
          inv[0][1] = - data_[0][1];
          inv[1][0] = - data_[1][0];
          inv[1][1] = data_[0][0];
          this->Determinant(det);
          inv *= 1.0 / det;
          break;

        case 3:
          // see Stoecker: "Taschenbuch Mathematischer Formeln und Moderner Verfahren" p.418
          inv.Resize(3,3);

          // === Old, recursive version ===
          // see Str: "Taschenbuch Mathematischer Formeln und Moderner Verfahren" p.418
          //for(UInt i=0; i<3; i++)
          //  for(UInt j=0; j<3; j++)
          //    inv[j][i] = Adjunct(i,j);
          this->Determinant(det);
          invDet = 1.0 / det;
          // === New, explicit version (from Wikipedia) ===
          inv[2][2] = (data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0])*invDet;
          inv[0][2] = (data_[0][1] * data_[1][2] - data_[0][2] * data_[1][1])*invDet;
          inv[1][2] = (data_[0][2] * data_[1][0] - data_[0][0] * data_[1][2])*invDet;

          inv[2][0] = (data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0])*invDet;
          inv[0][0] = (data_[1][1] * data_[2][2] - data_[1][2] * data_[2][1])*invDet;
          inv[1][0] = (data_[1][2] * data_[2][0] - data_[1][0] * data_[2][2])*invDet;

          inv[2][1] = (data_[0][1] * data_[2][0] - data_[0][0] * data_[2][1])*invDet;
          inv[1][1] = (data_[0][0] * data_[2][2] - data_[0][2] * data_[2][0])*invDet;
          inv[0][1] = (data_[0][2] * data_[2][1] - data_[0][1] * data_[2][2])*invDet;
          break;

        default:
          Double eps = 1e-40;
          TYPE pivot;
          TYPE pinv;

          //just copy the matrix
          inv.Resize(size_row_);
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
            inv.data_[0][k] = data_[0][k];

          //compute the inverse
          for ( UInt k=0; k<size_row_; k++) {
            pivot = inv[k][k];
            pinv  = 1.0 / pivot;
            // std::abs not possible for uint
            if ( std::abs(pivot) > eps ) {
              for (UInt j=0; j<size_row_; j++) {
                if ( j != k ) {
                  inv[k][j] = -inv[k][j] * pinv;
                  for ( UInt i=0; i<size_row_; i++ ) {
                    if ( i != k ) {
                      inv[i][j] = inv[i][j] + inv[i][k] * inv[k][j];
                    }
                  }
                }
              }
              for ( UInt i=0;  i<size_row_; i++ ) {
                inv[i][k] = inv[i][k] * pinv;
              }
              inv[k][k] = pinv;
              //	    std::cout << "inv current:\n" << inv << std::endl;
            }
            else {
              EXCEPTION("Get divison by zero in matrix inversion" );
            }
          }
        }
      }
    else {
      PseudoInvert(inv);
    }

  }

  template<> void Matrix<Integer>::Invert(Matrix<Integer> & inv) const {
    EXCEPTION("Integer matrix cannot be inverted");
  }

  template<> void Matrix<UInt>::Invert(Matrix<UInt> & inv) const {
    EXCEPTION("Integer matrix cannot be inverted");
  }

  template<class TYPE>
    void Matrix<TYPE>::Invert_Lapack() {
    EXCEPTION("General case not implemented");
  }
  
  template<> void Matrix<Complex>::Invert_Lapack() {
#ifdef CHECK_INDEX
    if( size_row_ != size_col_)
      EXCEPTION("Can only invert square matrices");
#endif

//TODO make sure this inversion is correct
    //std::cout<<"---------------------------------------------------\n"
    //		 <<"PLEASE TAKE CARE, THE INVERSION OF A COMPLEX MATRIX\n"
    //		 <<"USING LAPACK IS NOT THOROUGHLY TESTED!!!!!!!!!!!!!!\n"
    //		 <<"---------------------------------------------------"
	//		 <<std::endl;

    int *ipiv = new int[size_row_];
    int n = size_row_;
    int lwork = size_row_ * size_row_;
    std::complex<double> *work = new  std::complex<double>[lwork];
    int info;

    // calculate LU-factorization of block
    zgetrf(&n,&n,data_[0],&n,ipiv,&info);
    if( info != 0 ) {
      EXCEPTION("Error during LU-factorization of matrix. "
                << "Error value is " << info );
    }
    // invert matrix using previous LU factorization
    zgetri(&n,data_[0],&n,ipiv,work,&lwork,&info);
    if( info != 0 ) {
      EXCEPTION("Error during inversion of matrix. "
                << "Error value is " << info );
    }

    delete[] ipiv;
    delete[] work;
  }



  template<> void Matrix<Double>::Invert_Lapack() {
#ifdef CHECK_INDEX
    if( size_row_ != size_col_)
      EXCEPTION("Can only invert square matrices");
#endif

    int *ipiv = new int[size_row_];
    int n = size_row_;
    int lwork = size_row_ * size_row_;
    double *work = new double[lwork];
    int info;

    // calculate LU-factorization of block
    dgetrf(&n,&n,data_[0],&n,ipiv,&info);
    if( info != 0 ) {
      EXCEPTION("Error during LU-factorization of matrix. "
                << "Error value is " << info );
    }
    // invert matrix using previous LU factorization
    dgetri(&n,data_[0],&n,ipiv,work,&lwork,&info);
    if( info != 0 ) {
      EXCEPTION("Error during inversion of matrix. "
                << "Error value is " << info );
    }

    delete[] ipiv;
    delete[] work;
  }

  
  template<> void Matrix<Double>::Invert_Lapack(double & rcond, int & inf) {
#ifdef CHECK_INDEX
    if( size_row_ != size_col_) {
      EXCEPTION("Can only invert square matrices");
    }
#endif

    int *ipiv = new int[size_row_];
    int n = size_row_;
    int lwork = size_row_ * size_row_;
    double *work = new double[lwork];
    double *work1 = new double[4*n];
    int info1, info2, info3;

    // calculate LU-factorization of block
    dgetrf(&n,&n,data_[0],&n,ipiv,&info1);
    if( info1 != 0 ) {
      EXCEPTION("Error during LU-factorization of matrix. "
          << "Error value is " << info1 );
    }

    //compute 1-norm of the original matrix
    double anorm = 0.0;
    for(UInt i = 0; i < size_row_; ++i){
      for(UInt j = 0; j < size_col_; ++j){
        if(anorm < data_[i][j]) anorm = data_[i][j];
      }
    }

    char norm = '1'; //use the 1-norm, for infinity norm 'I'

    //compute the condition number
    dgecon(&norm, &n, data_[0], &n, &anorm, &rcond, work1, &n, &info2);
    if( info2 != 0 ) {
      EXCEPTION("Error during computation of condition number. "
          << "Error value is " << info2 );
    }

    // invert matrix using previous LU factorization
    dgetri(&n,data_[0],&n,ipiv,work,&lwork,&info3);
    if( info3 != 0 ) {
      EXCEPTION("Error during inversion of matrix. "
          << "Error value is " << info3 );
    }

    //check if any of the three operations above failed
    inf = 0;
    if(info1!=0 || info2!=0 || info3!=0) inf = 1;

    delete[] ipiv;
    delete[] work;
    delete[] work1;

  }

  template<class TYPE>
  bool Matrix<TYPE>::CholeskySolveLapack(Matrix<TYPE>& chol, Vector<TYPE>& x, const Vector<TYPE>& b, bool throw_exception)
  {
    throw Exception("CholeskySolveLapack is only implemented for double");
  }


  template<>
  bool Matrix<double>::CholeskySolveLapack(Matrix<double>& chol, Vector<double>& x, const Vector<double>& b, bool throw_exception)
  {
    char uplo = 'L';
    int n = size_row_;
    int nrhs = 1; // we could easily add more right hand sides!
    int info = 42;

    chol = *this; // will also resize and that stuff
    x = b; // need to copy content, as it will be overwritten

    // calculates Choleksy decomposition, stores lower part in chol and solves for x (overwritten
    dposv(&uplo, &n, &nrhs, chol.data_[0], &n, x.GetPointer(), &n, &info);

    // std::cout << "M:CSL b=" << b.ToString(3) << " x=" << x.ToString(3) << std::endl;
    // std::cout << "M:CSL info=" << info << " A=" << ToString(3) << std::endl;

    if(info == 0)
      return true;

    if(throw_exception)
      throw Exception("CholeskySolveLapack: dposv failed as matrix seems not to be positive definite, info=" + std::to_string(info));
    return false;
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::Adjunct (UInt i, UInt j) const
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("Undefined Matrix!" );
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_col_ ) 
      EXCEPTION("No quadratic matrix!" );
#endif
     
    assert(size_row_ == 3);
  
    Vector<Integer> iVec(2);
    Vector<Integer> jVec(2);
    UInt runningIndexI = 0;
    UInt runningIndexJ = 0;
  
    for (UInt actI = 0; actI<=2; actI++)
    {
      if (actI != i)
      {
        iVec[runningIndexI] = actI;
        runningIndexI++;
      }

      if (actI != j)
      {
        jVec[runningIndexJ] = actI;
        runningIndexJ++;
      }
    }

    // with pow(-1, i+j) a compiler complains: more than one instance of overloaded function "pow" matches the argument list  
    TYPE adj = (TYPE) std::pow(-1.0,(int) (i+j)) *
      ( data_[iVec[0]][jVec[0]] * data_[iVec[1]][jVec[1]] - 
        data_[iVec[0]][jVec[1]] * data_[iVec[1]][jVec[0]]);
    
    
    return adj;
    
  }
  
  template<class TYPE>
  Matrix<Double> Matrix<TYPE>::GetPart(  Global::ComplexPart part ) const {
    EXCEPTION( "Matrix::GetPart: Only Implemented for Real and Complex matrices!" );
    Matrix<Double> temp;
    return temp;
  }

  
  template<>
  Matrix<Double> Matrix<Double>::GetPart(  Global::ComplexPart part ) const {
    
    if ( part != Global::REAL ) {
      EXCEPTION("Matrix<Double>::GetPart: Only possible for REAL part." );
    }
      return *this;
    
  }

  template<>
  Matrix<Double> Matrix<Complex>::GetPart(Global::ComplexPart part) const
  {  
    Matrix<Double> ret(size_row_, size_col_);
    switch(part)
    {
    case Global::REAL:
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          ret[iRow][iCol]  = data_[iRow][iCol].real();
        }
      }
      break;
    case Global::IMAG:
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          ret[iRow][iCol]  = data_[iRow][iCol].imag();
        }
      }
      break;
    default:
      EXCEPTION("Matrix<Complex>::GetPart: Only possible for REAL or IMAG part!" );
    }
    
    return ret;
  }
  

  template<class TYPE>
  void Matrix<TYPE>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                              bool zeroOtherPart ) {
    EXCEPTION( "Matrix::SetPart: Only Implemented for Real and Complex matrices!" );
  }
  
  template<>
  void Matrix<Double>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                                bool zeroOtherPart)
  {
    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());
    assert((part == Global::REAL));

    *this = partMatrix;
  }

  template<>
  void Matrix<Complex>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                                 bool zeroOtherPart)
                                 {
    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());

    if( zeroOtherPart) {
      // ------------------
      //  Zero other part
      // ------------------
      switch(part)
      {
        case Global::REAL:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(partMatrix.data_[0][k], 0.0);
          }
          break;
        case Global::IMAG:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(0.0, partMatrix.data_[0][k]);
          }
          break;
        default:
          EXCEPTION( "Matrix<Complex>::SetPart: Only possible for REAL or IMAG part!" );
      }
    } else {


      // ------------------
      //  Keep other part
      // ------------------
      switch(part)
      {
        case Global::REAL:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(partMatrix.data_[0][k], data_[0][k].imag());
          }
          break;
        case Global::IMAG:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(data_[0][k].real(), partMatrix.data_[0][k]);
          }
          break;
        default:
          EXCEPTION( "Matrix<Complex>::SetPart: Only possible for REAL or IMAG part!" );
      }
    }
  }
  template<class TYPE>
  void Matrix<TYPE>::SetPartMult( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                                  Double factor, bool zeroOtherPart ) {
    EXCEPTION( "Matrix::SetPart: Only Implemented for Real and Complex matrices!" );
  }

  template<>
  void Matrix<Double>::SetPartMult( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                                    Double factor, bool zeroOtherPart)
  {
    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());
    assert((part == Global::REAL));

    *this = (partMatrix * factor);
  }

  template<> void Matrix<Complex>::
  SetPartMult( Global::ComplexPart part, const Matrix<Double> & partMatrix,
               Double factor, bool zeroOtherPart) {

    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());

    if( zeroOtherPart) {
      // ------------------
      //  Zero other part
      // ------------------
      switch(part)
      {
        case Global::REAL:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(partMatrix.data_[0][k] * factor, 0.0);
          }
          break;
        case Global::IMAG:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(0.0, partMatrix.data_[0][k] * factor);
          }
          break;
        default:
          EXCEPTION( "Matrix<Complex>::SetPart: Only possible for REAL or IMAG part!" );
      }
    } else {

      // ------------------
      //  Keep other part
      // ------------------
      switch(part)
      {
        case Global::REAL:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(partMatrix.data_[0][k] * factor, data_[0][k].imag());
          }
          break;
        case Global::IMAG:
          for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k) {
            data_[0][k] = Complex(data_[0][k].real(), partMatrix.data_[0][k] * factor);
          }
          break;
        default:
          EXCEPTION( "Matrix<Complex>::SetPart: Only possible for REAL or IMAG part!" );
      }
    }
  }

  // copies a submatrix at the position (row, col) into subMat, 
  // the amount of copied elements depends on the size of subMat
  template<class TYPE>
  void Matrix<TYPE>::GetSubMatrix(Matrix<TYPE>& subMat, UInt startRow, UInt startCol) const
  {

#ifdef CHECK_INITIALIZED
    if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
      EXCEPTION( "undefined matrix" );
#endif

#ifdef CHECK_INDEX
    if (((subMat.size_row_ + startRow) > size_row_) || ((subMat.size_col_ + startCol) > size_col_) )
      EXCEPTION( "Submatrx to be written is to large! " );
#endif

    for( UInt actRow=0; actRow < subMat.size_row_; actRow++)
      for( UInt actCol=0; actCol < subMat.size_col_; actCol++)
        subMat[actRow][actCol] = data_[actRow + startRow][actCol + startCol];  
  }



  // overwrites the matrix elements at the position (row, col) with subMat
  // in a rectangular (submatrix) way
  template<class TYPE>
  void Matrix<TYPE>::SetSubMatrix(const Matrix<TYPE>& subMat, UInt startRow, UInt startCol)
  {
#ifdef CHECK_INITIALIZED
    if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
      EXCEPTION( "undefined matrix" );
#endif
  
#ifdef CHECK_INDEX
    if ((subMat.size_row_ + startRow > size_row_) || (subMat.size_col_ + startCol > size_col_) )
      EXCEPTION("Submatrix to be read is to large!" );
#endif
  
    for( UInt actRow=0; actRow < subMat.size_row_; actRow++)
      for( UInt actCol=0; actCol < subMat.size_col_; actCol++)
        data_[actRow + startRow][actCol + startCol] = subMat[actRow][actCol];
  }



  // adds subMat to the matrix elements at the position (row, col)
  // in a rectangular (submatrix) way
  template<class TYPE>
  void Matrix<TYPE>::AddSubMatrix(const Matrix<TYPE>& subMat, UInt startRow, UInt startCol)
  {
#ifdef CHECK_INITIALIZED
    if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
      EXCEPTION("undefined matrix" );
#endif
  
#ifdef CHECK_INDEX
    if ((subMat.size_row_ + startRow > size_row_) || (subMat.size_col_ + startCol > size_col_) )
      EXCEPTION("Submatrix to be read is to large!");
#endif
  
    for( UInt actRow=0; actRow < subMat.size_row_; actRow++)
      for( UInt actCol=0; actCol < subMat.size_col_; actCol++)
        data_[actRow + startRow][actCol + startCol] += subMat[actRow][actCol];
  }



  /// converts a matrix into a vector, by appending successively all rows
  template<class TYPE>
  void Matrix<TYPE>::ConvertToVec_AppendRows(SingleVector & v) const
  {
    Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

    vec.Resize(size_row_ * size_col_);
  
    for( UInt i=0; i < size_row_; i++)
      for( UInt j=0; j < size_col_; j++)
        vec[i*(size_col_) + j] = data_[i][j];
  }

  /// converts a matrix into a vector, by appending successively all colums
  template<class TYPE>
  void Matrix<TYPE>::ConvertToVec_AppendCols(SingleVector &v) const
  {
    Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

    vec.Resize(size_row_ * size_col_);
  
    for( UInt actCol=0; actCol < size_col_; actCol++)
      for( UInt actRow=0; actRow < size_row_; actRow++)
          vec[actCol*size_row_ + actRow] = data_[actRow][actCol];
  }

  template<class TYPE>
  void Matrix<TYPE>::ConvertToVec_UpperTriangular(SingleVector & v) const
  {
    Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

    assert(size_row_ == size_col_);
    unsigned int size = size_row_;

    vec.Resize(((size * size - size)/2) + size);

    // 2*2 -> 11 22 12           = 00 11 01
    // 3*3 -> 11 22 33 23 13 12  = 00 11 22 12 02 01

    // diagonal first
    unsigned int pos = 0;
    for(unsigned int i = 0; i < size; i++)
      vec[pos++] = data_[i][i];

    // starting to write the upper triangular from behind upwards
    for(int c = size-1; c > 0; c--)
      for(int r = c-1; r >= 0; r--)
        vec[pos++] = data_[r][c];
  }


  /// scales the diagonal elements of a  matrix by a factor
  template<class TYPE>
  void Matrix<TYPE>::ScaleDiagElems(TYPE factor) 
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("Undefined Matrix!" );
#endif

#ifdef CHECK_INDEX 
    if (size_row_ != size_col_ ) 
      EXCEPTION("No square- matrix!" );
#endif

    for (UInt i = 0; i < size_row_; ++i)
      data_[i][i] *= factor;
  }


  template<class TYPE>
  void Matrix<TYPE>::GetRow(Vector<TYPE>& vec, UInt row) const
  {
    assert(size_row_ > row);
    vec.Resize(size_col_);

    for(unsigned int i = 0; i < size_col_; i++)
      vec[i] = (*this)[row][i]; // do it faster if you like
  }

  template<class TYPE>
  void Matrix<TYPE>::GetCol(Vector<TYPE>& vec, UInt col) const
  {
    assert(size_col_ > col);
    vec.Resize(size_row_);

    for(unsigned int i = 0; i < size_row_; i++)
      vec[i] = (*this)[i][col]; // do it faster if you like
  }

  template<class TYPE>
  void Matrix<TYPE>::SetCol(Vector<TYPE>& vec_in, UInt col) const
  {
      assert(size_col_ > col);
      for(unsigned int i = 0; i < size_row_; i++)
        data_[i][col] = vec_in[i];
  }

  template<class TYPE>
  void Matrix<TYPE>::GetColMin(Vector<TYPE>& vec) const
  {
    GetCol(vec, 0);

    for(unsigned int r = 0; r < size_row_; r++)
      for(unsigned int c = 1; c < size_col_; c++)
        vec[r] = std::min(vec[r], (*this)[r][c]); // do it faster if you like
  }

  template<>
  void Matrix<Complex>::GetColMin(Vector<Complex>& vec) const
  {
    GetCol(vec, 0);

    for(unsigned int r = 0; r < size_row_; r++)
      for(unsigned int c = 1; c < size_col_; c++) {
        vec[r].real(std::min(vec[r].real(), (*this)[r][c].real()));
        vec[r].imag(std::min(vec[r].imag(), (*this)[r][c].imag()));
      }
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::GetMax() const
  {
    Vector<TYPE> maxCol;
    GetColMax(maxCol);
    return maxCol.Max();
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::GetMin() const
  {
    Vector<TYPE> minCol;
    GetColMin(minCol);
    return minCol.Min();
  }

  template<class TYPE>
  void Matrix<TYPE>::GetColMax(Vector<TYPE>& vec) const
  {
    GetCol(vec, 0);

    for(unsigned int r = 0; r < size_row_; r++)
      for(unsigned int c = 1; c < size_col_; c++)
        vec[r] = std::max(vec[r], (*this)[r][c]); // do it faster if you like
  }

  template<>
  void Matrix<Complex>::GetColMax(Vector<Complex>& vec) const
  {
    GetCol(vec, 0);

    for(unsigned int r = 0; r < size_row_; r++)
      for(unsigned int c = 1; c < size_col_; c++) {
        vec[r].real(std::max(vec[r].real(), (*this)[r][c].real()));
        vec[r].imag(std::max(vec[r].imag(), (*this)[r][c].imag()));
      }
  }

  template<class TYPE>
  void Matrix<TYPE>::GetAbsValues(Matrix<TYPE>& AbsMatrix) const
  {
    assert(!(size_row_ == 0 || size_col_ == 0));

    AbsMatrix.Resize(size_row_, size_col_);

    for (UInt i = 0; i < size_row_; i++){
      for (UInt j = 0; j < size_col_; j++){
        AbsMatrix[i][j] = Abs(data_[i][j]);
      }
    }

  }

  /// gets the diagonal elements of a  matrix in a one column matrix
  template<class TYPE>
  void Matrix<TYPE>::GetDiagInMatrix(Matrix<TYPE>& columnMat) const
  {
    assert(!(size_row_ == 0 || size_col_ == 0));
    assert(size_row_ == size_col_ ); // check for square matrix

    columnMat.Resize(size_row_, 1);
    columnMat.Init();
  
    UInt i;
    for (i = 0; i < size_row_; i++)
      columnMat [i][0] = data_[i][i];

  }
  
  template<class TYPE>
  void Matrix<TYPE>::SetSubMatrixByInd( const Matrix<TYPE>& subMat,
                                        const StdVector<UInt>& rowInd,
                                        const StdVector<UInt>& colInd ) {

    UInt nrows = rowInd.GetSize();
    UInt ncols = colInd.GetSize();
    for( UInt iRow = 0; iRow < nrows; iRow++ ) {
      for( UInt iCol = 0; iCol < ncols; iCol++ ) {
#ifdef CHECK_INDEX
      if( rowInd[iRow] > size_row_ || 
          colInd[iCol] > size_col_ ) {
        EXCEPTION("Index pair (" << rowInd[iRow] << ", " << colInd[iCol]
                  << ") larger than matrix size being " << size_row_
                  << " x " << size_col_);
      }
#endif
      data_[rowInd[iRow]][colInd[iCol]]= subMat[iRow][iCol];
      }
    }
  }
  
  template<class TYPE>void Matrix<TYPE>::
  GetSubMatrixByInd( Matrix<TYPE>&ret, 
                     const StdVector<UInt>& rowInd,
                     const StdVector<UInt>& colInd ) const {

    UInt nrows = rowInd.GetSize();
    UInt ncols = colInd.GetSize();
    ret.Resize(nrows,ncols);
    for( UInt iRow = 0; iRow < nrows; iRow++ ) {
      for( UInt iCol = 0; iCol < ncols; iCol++ ) {
#ifdef CHECK_INDEX
      if( rowInd[iRow] > size_row_ || 
          colInd[iCol] > size_col_ ) {
        EXCEPTION("Index pair (" << rowInd[iRow] << ", " << colInd[iCol]
                  << ") larger than matrix size being " << size_row_
                  << " x " << size_col_);
      }
#endif
        ret[iRow][iCol] = data_[rowInd[iRow]][colInd[iCol]];
      }
    }
  }

  template<class TYPE>
  bool Matrix<TYPE>::IsSymmetric() const
  {
    if(!IsQuadratic())
      return false;

    for(UInt i = 1; i < size_row_; ++i)
      for(UInt j = i+1; j < size_col_; ++j)
        if(data_[i][j] != data_[j][i]) return false;
    
    return true;
  }

  // Putted into header
//
//  template<class TYPE>
//  bool Matrix<TYPE>::IsSymmetric(double eps) const
//  {
//    if(!IsQuadratic())
//      return false;
//
//    for(UInt i = 1; i < size_row_; ++i)
//      for(UInt j = i+1; j < size_col_; ++j)
//        if(!close(data_[i][j], data_[j][i]))
//          return false;
//
//    return true;
//  }

  template<class TYPE>
  TYPE Matrix<TYPE>::NormL2() const
  {
    TYPE result(0);
    
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      result += data_[0][k] * data_[0][k];

    return static_cast<TYPE>(std::sqrt(result)); // for compilers
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::Avg() const
  {
    assert(size_row_*size_col_ > 0);
    TYPE v(0);
    for(unsigned int i = 0, n = size_row_*size_col_; i < n; i++)
      v += data_[0][i];
    return v * (1./(size_row_*size_col_));
  }




  template<class TYPE>
  TYPE Matrix<TYPE>::DiffNormL2(const Matrix<TYPE>& other) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ != other.size_row_ || size_col_ != other.size_col_)
      EXCEPTION("Incompatible matrices");
#endif

    TYPE result(0);
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
    {
      TYPE tmp = data_[0][k] - other.data_[0][k];
      tmp *= tmp;
      result += tmp;
    }

    return static_cast<TYPE>(std::sqrt(result)); // for compilers
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::DiffNormL1(const Matrix<TYPE>& other) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ != other.size_row_ || size_col_ != other.size_col_)
      EXCEPTION("Incompatible matrices");
#endif

    TYPE result(0);
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      result += static_cast<TYPE>(Abs(data_[0][k] - other.data_[0][k]));

    return result;
  }

  template<class TYPE>
  bool Matrix<TYPE>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if(std::isnan(data_[0][k])) return true;

    return false;
  }

  template<>
  bool Matrix<int>::ContainsNaN() const
  {
    return false;
  }

  template<>
  bool Matrix<unsigned int>::ContainsNaN() const
  {
    return false;
  }


  template<>
  bool Matrix<Complex>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
    {
      if(std::isnan(data_[0][k].real())) return true;
      if(std::isnan(data_[0][k].imag())) return true;
    }
    return false;
  }

  template<class TYPE>
  bool Matrix<TYPE>::ContainsInf() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if(std::isinf(data_[0][k])) return true;

    return false;
  }

  template<>
  bool Matrix<int>::ContainsInf() const
  {
    return false;
  }

  template<>
  bool Matrix<unsigned int>::ContainsInf() const
  {
    return false;
  }

  template<>
  bool Matrix<Complex>::ContainsInf() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
    {
      if(std::isinf(data_[0][k].real())) return true;
      if(std::isinf(data_[0][k].imag())) return true;
    }
    return false;
  }


  // Alternate version of symmetry checker, that will report
  // asymmetries to standard output

  // template<class TYPE> bool Matrix<TYPE>::IsSymmetric() {
  //   bool amSymm = true;
  //   for ( UInt i = 1; i < size_row_; i++ ) {
  //     for ( UInt j = i+1; j < size_col_; j++ ) {
  //       if ( data_[i][j] != data_[j][i] ) {
  //      amSymm = false;
  //      (*debug) << " (" << i << "," << j << "): " << data_[i][j] << " <--> "
  //               << data_[j][i] << "(abs.diff= " << (data_[i][j]-data_[j][i])
  //               << ") "
  //               << "(rel.diff= "
  //               << (data_[i][j]-data_[j][i])/data_[i][j] << ")\n";
  //       }
  //     }
  //   }
  //   return amSymm;
  // }

  template<class TYPE>
  Matrix<TYPE> Conj( const Matrix<TYPE>& m) {
    return m;
  }

  template<>
  Matrix<Complex> Conj( const Matrix<Complex>& m) {
    UInt numRows = m.GetNumRows();
    UInt numCols = m.GetNumCols();
    Matrix<Complex> ret(numRows, numCols);
    for( UInt i = 0; i < numRows; i++ ) {
      for (UInt j = 0; j < numCols; j++ ) {
        ret[i][j] = std::conj(m[i][j]);
      }
    }
    return ret;
  }

  template<class TYPE>
  Matrix<TYPE> Herm( const Matrix<TYPE>& m ) {
    return m;
  }

  template<>
  Matrix<Complex> Herm( const Matrix<Complex>& m ) {
    UInt numRows = m.GetNumRows();
    UInt numCols = m.GetNumCols();
    Matrix<Complex> herm(numCols, numRows);
    for( UInt i = 0; i < numCols; i++ ) {
      for (UInt j = 0; j < numRows; j++ ) {
        herm[i][j] = std::conj(m[j][i]);
      }
    }
    return herm;
  }

  template Matrix<Double> Conj( const Matrix<Double>& m);
  template Matrix<Double> Herm( const Matrix<Double>& m);

#ifdef __GNUC__
  template class Matrix<Double>;
  template class Matrix<Integer>;
  template class Matrix<UInt>;
  template class Matrix<Complex>;
#endif

#if defined (_MSC_VER) 
  template class Matrix<Double>;
  template class Matrix<Integer>;
  template class Matrix<UInt>;
  template class Matrix<Complex>;
  template void Matrix<Double>::Add(Double, const Matrix<Double>&);
  template void Matrix<Integer>::Add(Integer, const Matrix<Integer>&);
  template void Matrix<UInt>::Add(UInt, const Matrix<UInt>&);
  template void Matrix<Complex>::Add(Complex, const Matrix<Complex>&);
//  template void Matrix<Complex>::Add(Complex, const Matrix<Double>&);
  template void Matrix<Double>::AddT(Double, const Matrix<Double>&);
  template void Matrix<Integer>::AddT(Integer, const Matrix<Integer>&);
  template void Matrix<UInt>::AddT(UInt, const Matrix<UInt>&);
  template void Matrix<Complex>::AddT(Complex, const Matrix<Complex>&);
//  template void Matrix<Complex>::AddT(Complex, const Matrix<Double>&);
  /*
*/
#endif


// explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Matrix<UInt>
#pragma instantiate Matrix<Integer>
#pragma instantiate Matrix<Double>
#pragma instantiate Matrix<Complex>
#endif

} // end of namespace

