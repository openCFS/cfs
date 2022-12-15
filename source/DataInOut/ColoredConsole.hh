//------------------------------------------------------------------------------
// Console.h: interface for the Console manipulators.
//------------------------------------------------------------------------------

#if !defined( CONSOLE_MANIP_H__INCLUDED )
#define CONSOLE_MANIP_H__INCLUDED

//------------------------------------------------------------------------------

//------------------------------------------------------------------"includes"--
#include <iostream>
#include <iomanip>

#if defined(WIN32)
#define NOMINMAX
#include <windows.h>
#endif

namespace CoupledField
{

#if defined(WIN32)
  static const unsigned int bgMask( BACKGROUND_BLUE      |
                            BACKGROUND_GREEN     |
                            BACKGROUND_RED       |
                            BACKGROUND_INTENSITY   );
  static const unsigned int fgMask( FOREGROUND_BLUE      |
                            FOREGROUND_GREEN     |
                            FOREGROUND_RED       |
                            FOREGROUND_INTENSITY   );

  static const unsigned int fgBlack    ( 0 );
  static const unsigned int fgLoRed    ( FOREGROUND_RED   );
  static const unsigned int fgLoGreen  ( FOREGROUND_GREEN );
  static const unsigned int fgLoBlue   ( FOREGROUND_BLUE  );
  static const unsigned int fgLoCyan   ( fgLoGreen   | fgLoBlue );
  static const unsigned int fgLoMagenta( fgLoRed     | fgLoBlue );
  static const unsigned int fgLoYellow ( fgLoRed     | fgLoGreen );
  static const unsigned int fgLoWhite  ( fgLoRed     | fgLoGreen | fgLoBlue );
  static const unsigned int fgGray     ( fgBlack     | FOREGROUND_INTENSITY );
  static const unsigned int fgHiWhite  ( fgLoWhite   | FOREGROUND_INTENSITY );
  static const unsigned int fgHiBlue   ( fgLoBlue    | FOREGROUND_INTENSITY );
  static const unsigned int fgHiGreen  ( fgLoGreen   | FOREGROUND_INTENSITY );
  static const unsigned int fgHiRed    ( fgLoRed     | FOREGROUND_INTENSITY );
  static const unsigned int fgHiCyan   ( fgLoCyan    | FOREGROUND_INTENSITY );
  static const unsigned int fgHiMagenta( fgLoMagenta | FOREGROUND_INTENSITY );
  static const unsigned int fgHiYellow ( fgLoYellow  | FOREGROUND_INTENSITY );
  static const unsigned int bgBlack    ( 0 );
  static const unsigned int bgLoRed    ( BACKGROUND_RED   );
  static const unsigned int bgLoGreen  ( BACKGROUND_GREEN );
  static const unsigned int bgLoBlue   ( BACKGROUND_BLUE  );
  static const unsigned int bgLoCyan   ( bgLoGreen   | bgLoBlue );
  static const unsigned int bgLoMagenta( bgLoRed     | bgLoBlue );
  static const unsigned int bgLoYellow ( bgLoRed     | bgLoGreen );
  static const unsigned int bgLoWhite  ( bgLoRed     | bgLoGreen | bgLoBlue );
  static const unsigned int bgGray     ( bgBlack     | BACKGROUND_INTENSITY );
  static const unsigned int bgHiWhite  ( bgLoWhite   | BACKGROUND_INTENSITY );
  static const unsigned int bgHiBlue   ( bgLoBlue    | BACKGROUND_INTENSITY );
  static const unsigned int bgHiGreen  ( bgLoGreen   | BACKGROUND_INTENSITY );
  static const unsigned int bgHiRed    ( bgLoRed     | BACKGROUND_INTENSITY );
  static const unsigned int bgHiCyan   ( bgLoCyan    | BACKGROUND_INTENSITY );
  static const unsigned int bgHiMagenta( bgLoMagenta | BACKGROUND_INTENSITY );
  static const unsigned int bgHiYellow ( bgLoYellow  | BACKGROUND_INTENSITY );

  static const unsigned int fgbgReset ( (bgMask  | fgMask) + 1);

#endif

  static class ColoredConsole
  {
#if defined(WIN32)
  private:
    HANDLE                      hCon;
    unsigned int                        cCharsWritten;
    CONSOLE_SCREEN_BUFFER_INFO  csbi;
    unsigned int                        dwConSize;
	WORD                        wDefaultAttributes;
#endif

  public:
    static bool colorise;
    static bool suppressed;

  public:
    ColoredConsole() {
#if defined(WIN32)
      hCon = GetStdHandle( STD_OUTPUT_HANDLE );
      GetConsoleScreenBufferInfo( hCon, &csbi );
	  wDefaultAttributes = csbi.wAttributes;
#endif
    }

  private:
    void GetInfo()
    {
#if defined(WIN32)
      GetConsoleScreenBufferInfo( hCon, &csbi );
      dwConSize = csbi.dwSize.X * csbi.dwSize.Y;
   #endif
    }
  public:
    void Clear()
    {
#if defined(WIN32)
      COORD coordScreen = { 0, 0 };

      GetInfo();
      FillConsoleOutputCharacter( hCon, TEXT(' '),
                                  dwConSize,
                                  coordScreen,
                                  (DWORD*) &cCharsWritten );
      GetInfo();
      FillConsoleOutputAttribute( hCon,
                                  csbi.wAttributes,
                                  dwConSize,
                                  coordScreen,
                                  (DWORD*) &cCharsWritten );
      SetConsoleCursorPosition( hCon, coordScreen );
   #endif
    }

    void SetColor( unsigned int wRGBI,
                   unsigned int Mask,
                   std::string modif,
                   std::ostream& os )
    {
      if(suppressed)
        return;

      if(!colorise)
        return;

#if defined(WIN32)
      if(&os == &std::cout)
        hCon = GetStdHandle( STD_OUTPUT_HANDLE );
      else if(&os == &std::cerr)
        hCon = GetStdHandle( STD_ERROR_HANDLE );
      else
        return;

      GetInfo();
      if((wRGBI == fgbgReset) && (Mask == fgbgReset))
	  {
        csbi.wAttributes = wDefaultAttributes;
	  }
	  else
	  {
        csbi.wAttributes &= WORD(Mask);
        csbi.wAttributes |= WORD(wRGBI);
	  }

	  SetConsoleTextAttribute( hCon, csbi.wAttributes );
   #else
      os << modif;
   #endif

    }
  } console;

  inline std::ostream& clr( std::ostream& os )
  {
    os.flush();
    console.Clear();
    return os;
  }

  inline std::ostream& fg_red( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgHiRed, bgMask, "", os );
 #else
    // light red    
    // console.SetColor( 0, 0, "\033[1;31m", os);
    console.SetColor( 0, 0, "\033[01;31m", os);
 #endif
    return os;
  }

  inline std::ostream& fg_green( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgHiGreen, bgMask, "", os  );
 #else
    // light green    
    // console.SetColor( 0, 0, "\033[32;1m", os);
    console.SetColor( 0, 0, "\033[01;32m", os);
 #endif

    return os;
  }

  inline std::ostream& fg_blue( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgHiBlue, bgMask, "", os  );
 #else
    // light blue   
    // console.SetColor( 0, 0, "\033[34;1m", os);
    console.SetColor( 0, 0, "\033[01;34m", os);
 #endif

    return os;
  }

  inline std::ostream& fg_white( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgHiWhite, bgMask, "", os  );
 #else
    console.SetColor( 0, 0, "\033[37;1m", os);
 #endif

    return os;
  }

  inline std::ostream& fg_cyan( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgHiCyan, bgMask, "", os  );
 #else
    console.SetColor( 0, 0, "\033[36;1m", os);
 #endif

    return os;
  }

  inline std::ostream& fg_magenta( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgHiMagenta, bgMask, "", os );
 #else
    console.SetColor( 0, 0, "\033[35;1m", os);
 #endif

    return os;
  }

  inline std::ostream& fg_yellow( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgHiYellow, bgMask, "", os  );
 #else
    console.SetColor( 0, 0, "\033[33;1m", os);
 #endif

    return os;
  }

  inline std::ostream& fg_black( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgBlack, bgMask, "", os  );
 #else
    console.SetColor( 0, 0, "\033[0m", os);
 #endif

    return os;
  }

  inline std::ostream& fg_gray( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgGray, bgMask, "", os  );
 #else
    console.SetColor( 0, 0, "\033[30;1m", os);
 #endif

    return os;
  }

  inline std::ostream& bg_red( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgHiRed, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_green( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgHiGreen, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_blue( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgHiBlue, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_white( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgHiWhite, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_cyan( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgHiCyan, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_magenta( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgHiMagenta, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_yellow( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgHiYellow, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_black( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgBlack, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& bg_gray( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( bgGray, fgMask, "", os  );
 #else
 #endif

    return os;
  }

  inline std::ostream& fg_reset( std::ostream& os )
  {
    os.flush();
#if defined(WIN32)
    console.SetColor( fgbgReset, fgbgReset, "", os  );
 #else
    console.SetColor( 0, 0, "\033[0m", os);
 #endif

    return os;
  }

}

 // HACK to prevent name clashes while compiling WriteInfo.cc
#undef CONST
#undef CreateFile

//------------------------------------------------------------------------------
#endif //!defined ( CONSOLE_MANIP_H__INCLUDED )

// i:\dev\Intel\Compiler\C++\10.0.025\Ia32\Bin\icl.exe /E /nologo /TP -DWIN32 /W3 /Zm1000 /GX /GR /debug:all /W1 /Wcheck /Qdiag-disable:654,1125     /DBOOST_ALL_NO_LIB /DBOOST_ALL_DYN_LINK /MDd /Zi /Od /GZ -Ii:\dev\CFSDEPS\include -Ii:\NACS\CACHE -Ii:\NACS\PATH -Ii:\NACS\INTERNAL -Ii:\NACS\CFSDEPS_INCLUDE_DIR -Id:\NACS_BUILD_INTEL_NATIVE\include -Ii:\NACS\source -Ii:\NACS\source\OLAS -c i:\NACS\source\DataInOut\WriteInfo.cc
