// Definition of the GSISocket class

#ifndef GSISocket_class
#define GSISocket_class


#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <string>
#include <arpa/inet.h>

const int MAXHOSTNAME = 200;
const int MAXCONNECTIONS = 1;
const int MAXRECV = 500;

class GSISocket
{
 public:
  GSISocket(int timeout=-1);
  virtual ~GSISocket();

  // Server initialization
  bool create();
  bool bind ( const int port );
  bool listen() const;
  int accept ( GSISocket& );

  // Client initialization
  bool connect ( const std::string& host, const int port );

  FILE* getReadHandle();
  FILE* getWriteHandle();

  void set_non_blocking ( const bool );

  bool is_valid() const { return m_sock != -1; }

  int read(void* buf, int nbytes, int timeout);
  int write(void* buf, int nbytes, int timeout);
  
 private:

  int m_sock;
  sockaddr_in m_addr;
  FILE* file_read;
  FILE* file_write;

 protected:
  void cleanup();
  void shutdownServer();
  int mayDoIO(int type, int timeout);
  FILE* sock_to_file(int sock, const char *mode);

  int timeout_;
};


#endif
