#ifdef USE_RCSID
static const char RCSid_GSISocket[] = "$Id$";
#endif

/*----------------------------------------------------------------------
  |
  |
  | $Log$
  | Revision 1.4  2005/05/23 22:15:31  ahauck
  | Corrected indentation of all cc/hh-files according to xemacs' 'gnu' style
  | and according to our coding rules, so all TABs were replaced by spaces.
  |
  | Revision 1.3  2005/05/18 19:26:03  strieben
  | Upgraded GSI library to newest available version.
  |
  | Revision 1.1.1.1  2004/08/31 15:53:00  simon
  | Initial GSI import
  |
  |
  +---------------------------------------------------------------------*/

#include <cstring>
#include <iostream>

#include "GSISocket.hh"

namespace GridlibSocketInterface
{

  static void file_status(int32 fd) {
    int32 open_modus, wert;

    if((wert=fcntl(fd, F_GETFL, 0)) == -1) {
      printf("Fehler bei fcntl\n");
    }

    open_modus = wert & O_ACCMODE;

    if(open_modus == O_RDONLY) printf("read only");
    else if(open_modus == O_WRONLY) printf("write only");
    else if(open_modus == O_RDWR) printf("read write");
    else printf("unbekannter open modus");

    if(wert & O_APPEND) printf(", append");
    if(wert & O_NONBLOCK) printf(", nonblocking");
#ifdef O_SYNC
    if(wert & O_SYNC) printf(", O_SYNC gesetzt");
#endif
    printf("\n");
  }


  Socket::Socket(int32 timeout) :
    m_sock ( -1 ), file_read(NULL), file_write (NULL), timeout_(timeout)
  {
    //  file_read = NULL;
    //  file_write = NULL;

    memset ( &m_addr,
             0,
             sizeof ( m_addr ) );

  }

  Socket::~Socket()
  {
    cleanup();
  }

  bool Socket::create()
  {
    m_sock = socket ( AF_INET,
                      SOCK_STREAM,
                      0 );

    if ( ! is_valid() )
      return false;


    // TIME_WAIT - argh
    int32 on = 1;
    if ( setsockopt ( m_sock, SOL_SOCKET, SO_REUSEADDR, ( const int8* ) &on, sizeof ( on ) ) == -1 )
      return false;


    return true;
  }



  bool Socket::bind ( const int32 port )
  {

    if ( ! is_valid() )
      {
        return false;
      }
  


    m_addr.sin_family = AF_INET;
    m_addr.sin_addr.s_addr = INADDR_ANY;
    m_addr.sin_port = htons ( port );

    int32 bind_return = ::bind ( m_sock,
                                 ( struct sockaddr * ) &m_addr,
                                 sizeof ( m_addr ) );
  

    if ( bind_return == -1 )
      {
        return false;
      }
  
    return true;
  }


  bool Socket::listen() const
  {
    if ( ! is_valid() )
      {
        return false;
      }

    int32 listen_return = ::listen ( m_sock, MAXCONNECTIONS );


    if ( listen_return == -1 )
      {
        return false;
      }

    return true;
  }


  int32 Socket::accept ( Socket& new_socket )
  {
    int32 addr_length = sizeof ( m_addr );
    int32 ret;

    fd_set acceptset;
    struct timeval tv;
  
    if(timeout_ >= 0) 
      {
      
        tv.tv_sec = timeout_;
        tv.tv_usec = 0;

        FD_ZERO(&acceptset);
      
        FD_SET(m_sock, &acceptset);
  
        if((ret = select(m_sock+1, &acceptset, NULL, NULL, &tv)) < 0) 
          {
            if(errno == EINTR)
              return 0;
            else
              return -1;
          }
        else 
          {
            if(ret == 0) 
              {
                //              std::cout << "No activity on accept socket for 5 seconds" << std::endl;
                return 0;
              }
          }
      }

    // Here we can be sure that someone wants to connect
    new_socket.m_sock = ::accept ( m_sock, ( sockaddr * ) &m_addr, (socklen_t*) &addr_length );
    if ( new_socket.m_sock <= 0 ) 
      return -2;
    else {
      //      printf("Connection accepted\n");
      //              file_status(new_socket.m_sock);
      
      if((new_socket.file_read = sock_to_file(new_socket.m_sock, "rb")) == NULL)
        //  fdopen( new_socket.m_sock, "rb" ) ) == NULL ) 
        {
          new_socket.cleanup();
          
          return -3;
        }
      if( ( new_socket.file_write = sock_to_file(dup(new_socket.m_sock), "wb")) == NULL)
        // fdopen( dup(new_socket.m_sock), "wb" ) ) == NULL ) 
        {
          new_socket.cleanup();
          
          return -4;
        }
      
      /*
        if( setvbuf(new_socket.file_read, (int8 *)NULL, _IONBF, 0) != 0)
        {
        new_socket.cleanup();
          
        return -5;
        }

        if( setvbuf(new_socket.file_write, (int8 *)NULL, _IONBF, 0) != 0)
        {
        new_socket.cleanup();
          
        return -6;
        }
      */

      return 1;
    }
  }

  int32 Socket::write(void* buf, int32 nbytes, int32 timeout)
  {
    int32 ret = -1;

    if(!is_valid())    
      return -1;  

    // we want to write data to the socket
    ret = mayDoIO(0, timeout);
  
    //  std::cout << "Socket write ret1: " << ret << std::endl;
  
    if(ret <= 0)
      return ret;

    register const int8 *ptr = (const int8 *) buf;
    register int32 bytesleft = nbytes;

    do
      {
        register int32 rc;
      
        do
          rc = ::write(m_sock, ptr, bytesleft);
        while (rc < 0 && errno == EINTR);
      
        if (rc < 0)
          return -1;
      
        bytesleft -= rc;
        ptr += rc;
      }
    while (bytesleft);
  
    return nbytes;



    /*
      socklen_t sndbufsize, optlen;

      getsockopt(m_sock, SOL_SOCKET, SO_SNDBUF, &sndbufsize, &optlen);
  
      std::cout << "Socket sndbufsize: " << sndbufsize << std::endl;

      // we want to write data to the socket
      ret = mayDoIO(0, timeout);
  
      std::cout << "Socket write ret1: " << ret << std::endl;
  
      if(ret <= 0)
      return ret;
  
      ret = ::write(m_sock, buf, nbytes);
  
      perror("CFSSocket");
  
      std::cout << "Socket write ret2: " << ret << std::endl;
      if(ret != nbytes)
      return -1;
      else
      return ret;
    */
  }

  int32 Socket::read(void* buf, int32 nbytes, int32 timeout)
  {
    int32 ret = -1;
  
    if(!is_valid())    
      return -1;
  
    // we want to read data from the socket
    ret = mayDoIO(1, timeout);
  
    if(ret <= 0)
      return ret;
  
    register int8 *ptr = (int8 *) buf;
    register int32 bytesleft = nbytes;
  
    do
      {
        register int32 rc;
      
        do
          rc = ::read(m_sock, ptr, bytesleft);
        while (rc < 0 && errno == EINTR);
      
        if (rc == 0)
          return nbytes - bytesleft;
        else if (rc < 0)
          return -1;
      
        bytesleft -= rc;
        ptr += rc;
      }
    while (bytesleft);
  
    return nbytes;


    /*
      ret = ::read(m_sock, buf, nbytes);

      if(ret != nbytes)
      return -1;
      else
      return ret;
    */
  }

  void Socket::cleanup() 
  {
    //  std::cout << "cleaning up file_read" << std::endl;
    if(file_read)
      while(fclose(file_read) == -1) 
        {
          if(errno == EINTR)
            continue;
        }
    file_read = NULL;
  
    //  std::cout << "cleaning up file_write" << std::endl;
    if(file_write)
      while(fclose(file_write) == -1) 
        {
          if(errno == EINTR)
            continue;
        }
    file_write = NULL;
  
    // according to the man page of fdopen the following should not be necessary
    //  if ( is_valid() )
    //    ::close ( m_sock );
  
    if ( is_valid() )
      ::close ( m_sock );
  
    m_sock = -1;
  

    memset ( &m_addr,
             0,
             sizeof ( m_addr ) );
  }

  void Socket::shutdownServer() 
  {
    if(is_valid())
      while(close(m_sock) == -1) 
        {
          if(errno == EINTR)
            continue;
        }
  
    m_sock = -1;
  
  
    memset ( &m_addr,
             0,
             sizeof ( m_addr ) );
  }

  FILE* Socket::getReadHandle() {
    //  if(file_read==NULL) printf("WARNING: NULL readhandle\n");
    return file_read;
  }

  FILE* Socket::getWriteHandle() {
    //  if(file_read==NULL) printf("WARNING: NULL writehandle\n");
    return file_write;
  }

  bool Socket::connect ( const std::string& hostname, const int32 port )
  {
    if ( ! is_valid() ) return false;

    struct in_addr inadr;
    struct hostent *host;
  
    if(inet_aton(hostname.c_str(), &inadr)) 
      host = gethostbyaddr((int8*) &inadr, sizeof(inadr), AF_INET);
    else
      host = gethostbyname(hostname.c_str());
  
    if(!host) return false;
  

    m_addr.sin_family = AF_INET;
    m_addr.sin_port = htons ( port );
    memcpy(&m_addr.sin_addr, host->h_addr_list[0], sizeof(m_addr.sin_addr));
  
    //  int status = inet_pton ( AF_INET, host.c_str(), &m_addr.sin_addr );
    //  if ( errno == EAFNOSUPPORT ) return false;

    int32 status;

    status = ::connect ( m_sock, ( sockaddr * ) &m_addr, sizeof ( m_addr ) );

    if ( status == 0 ) {
      //    file_status(m_sock);
      if( ( file_read = sock_to_file( m_sock, "rb" ) ) == NULL ) {
        cleanup();
        return false;
      }
      if( ( file_write = sock_to_file( dup(m_sock), "wb" ) ) == NULL ) {
        cleanup();
        return false;
      }

      return true;
    }
    else
      return false;
  }

  void Socket::set_non_blocking ( const bool b )
  {

    int32 opts;

    opts = fcntl ( m_sock,
                   F_GETFL );

    if ( opts < 0 )
      {
        return;
      }

    if ( b )
      opts = ( opts | O_NONBLOCK );
    else
      opts = ( opts & ~O_NONBLOCK );

    fcntl ( m_sock,
            F_SETFL,opts );

  }

  int32 Socket::mayDoIO(int32 type, int32 timeout)
  {
    int32 ret = -1;

    fd_set fdset;
    struct timeval tv;

    if(!is_valid())
      return -1;

    if(timeout < 0)
      return 1;
    else 
      {
        tv.tv_sec = timeout;
        tv.tv_usec = 0;
      }
  
    FD_ZERO(&fdset);
  
    FD_SET(m_sock, &fdset);
  
    //  std::cout << "Write descriptor " << desc_write << std::endl;

    if(type == 0)
      ret = select(m_sock+1, NULL, &fdset, NULL, &tv);
    if(type == 1)
      ret = select(m_sock+1, &fdset, NULL, NULL, &tv);

    if(ret < 0) 
      {
        if(errno == EINTR)
          {
            return 0;
          }
        else 
          {
            return -1;
          }
      }
    else 
      {
        if(ret == 0) 
          {
            return 0;
          }
        else
          {
            return 1;
          }
      }
  }


#define SOCK_SMALLBUF (4*1024)
#define SOCK_LARGEBUF (28*1024)

  FILE *Socket::sock_to_file(int32 sock, const int8 *mode)
  {
    int32 bufsz = SOCK_SMALLBUF;
    int32 bufmode = _IOFBF;

    FILE *fp = fdopen(sock, (*mode == 'w') ? "w" : "r");

    switch (mode[1])
      {
      case 't':
        bufmode = _IOLBF;
        break;

      case 'b':
        bufsz = SOCK_LARGEBUF * 2;
        setsockopt(sock, SOL_SOCKET,
                   (*mode == 'w') ? SO_SNDBUF : SO_RCVBUF,
                   &bufsz, sizeof(bufsz));
        bufsz = SOCK_LARGEBUF;
        break;
      }

    setvbuf(fp, NULL, bufmode, bufsz);

    return fp;
  }
 
}
