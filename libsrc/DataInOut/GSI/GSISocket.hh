/*----------------------------------------------------------------------
  |
  | $Id$
  |
  +---------------------------------------------------------------------*/


#ifndef GSI_SOCKET
#define GSI_SOCKET

/* According to earlier standards */
#include <sys/time.h>
#include <sys/select.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <string>

#include "GSITypeDefs.hh"

namespace GridlibSocketInterface
{

  const int32 MAXHOSTNAME = 200;
  const int32 MAXCONNECTIONS = 1;
  const int32 MAXRECV = 500;

  class Socket
  {
  public:
    Socket(int32 timeout=-1);
    virtual ~Socket();

    // Server initialization
    bool create();
    bool bind ( const int32 port );
    bool listen() const;
    int32 accept ( Socket& );

    // Client initialization
    bool connect ( const std::string& host, const int32 port );

    FILE* getReadHandle();
    FILE* getWriteHandle();

    void set_non_blocking ( const bool );

    bool is_valid() const { return m_sock != -1; }

    int32 read(void* buf, int32 nbytes, int32 timeout);
    int32 write(void* buf, int32 nbytes, int32 timeout);

  private:

    int32 m_sock;
    sockaddr_in m_addr;
    FILE* file_read;
    FILE* file_write;

  protected:
    void cleanup();
    void shutdownServer();
    int32 mayDoIO(int32 type, int32 timeout);
    FILE* sock_to_file(int32 sock, const int8 *mode);

    int32 timeout_;
  };
 
}

#endif //GSI_SOCKET
