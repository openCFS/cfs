#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

#define DEFAULT_BITRATE 5000000.0

int xres,yres,start,stop,add;
double maxMB;

char S[5000];

// this function is auxiliary for function get
std::string::size_type getpos(std::ifstream & in, const std::string keyword)
{
  std::string::size_type help,pos=std::string::npos;
  std::string buf;

  in.seekg(0, std::ios::beg);
   while ( pos == std::string::npos && !in.eof() )
  {
    help=in.tellg();
    std::getline(in, buf, '\n');
    pos=buf.find(keyword);
  }

  if (in.eof()) std::cerr << "Can't find keyword " << keyword  << std::endl;

  pos=buf.find("=");

  return pos+help+1;
}

// we use this function for reading parameters from mkmpeg.conf file
template<class TypeVal>
void get(std::ifstream & ainfile, const std::string keyword, TypeVal & val)
{
   std::string::size_type pos,pos1=0;

   pos=getpos(ainfile,keyword);

   ainfile.seekg(pos, std::ios::beg);
   ainfile >> val;
}

int main(int argc,char **argv)
{
  // introduction: usage of programm 
   if(argc!= 4 && argc!=5)
     {
       std::cout << "\033[31mUsage\033[0m: mkmpeg name start stop <add>" << std::endl
              << "\t \033[36m name \033[0m: name of input sequence of gmv-files" << std::endl
       << "\t \033[36m start \033[0m: number of first file " << std::endl
       << "\t \033[36m stop \033[0m: number of last file" << std::endl << std::endl
       << "\t \033[36m <add> \033[0m: step size between data files, by default - 1" << std::endl << std::endl;
        return 1;
     }

// open configuration file
   std::ifstream infile;
   infile.open("gmv2mpeg.config");
   if (!infile.good()) { std::cerr << "Error: can't find file gmv2mpeg.config" << std::endl; exit(1);}

   std::string attr;

// read parameters
   get(infile,"xresolution",xres);
   get(infile,"yresolution",yres);
   get(infile,"max_MB",maxMB);
   get(infile,"attrfile",attr);

//check that attribute file exits
   sprintf(S,"test -r %s",attr.c_str());
   if(system(S))
          {
             printf("Can't find file %s\n",attr.c_str());
             exit(1);
          }

   start=atoi(argv[2]);
   stop=atoi(argv[3]);

   if(argc==5) add=atoi(argv[4]);
        else add=1;

// create snapshots in format YUV
   int i,counter;
   for(i=start;i<=stop;i+=add)
     {
        sprintf(S,"test -r %s.gmv00%i",argv[1],i);
        if(system(S))
          {
             printf("Can't find file  %s.gmv00%i!\n",argv[1],i);
             exit(1);
          }
        sprintf(S,"gmv -m -a %s -w 0 0 %d %d -i %s.gmv00%d -s\n", attr.c_str(),xres,yres, argv[1],i);
        printf(S);
        system(S);
        sprintf(S,"sgitopnm AzsnapgmvAz | ppmtoyuvsplit %s%d\n",argv[1],i+1);
        printf(S);
        system(S);
        printf("Frame %i is done!\n",i+1);
        counter++;
     }
   remove("AzsnapgmvAz");

// define number of bits
   long int max_bits;
   if(maxMB==0) max_bits=(long int)(DEFAULT_BITRATE*counter)/25;
    else max_bits=(long int)maxMB*1024*1024*8;

   std::cout << " Start mpeg encoding " << std::endl;

// do mpeg encoding
   sprintf(S,"mpeg -PF -p 3 -a 1 -b %d -h %d -v %d -x %ld -s %s.mpeg %s >/dev/null",counter,xres,yres,max_bits,argv[1],argv[1]);

   puts(S);
   system(S);

   std::string savesnapshot;
   get(infile, "save_snapshot",savesnapshot);
   if (savesnapshot == "no")
   {
   sprintf(S,"rm -f %s*U; rm -f %s*V; rm -f %s*Y", argv[1], argv[1], argv[1]);
   system(S);
   }

   return 1;
}
