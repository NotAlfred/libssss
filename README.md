# libssss
This is a port of [point-at-infinity](http://point-at-infinity.org/ssss/index.html)'s ssss implementation turned into a reusable C library
## Build
### Ubuntu
```bash
sudo apt install libgmp-dev
git clone https://github.com/NotAlfred/libssss.git
cd libssss
make
```

## Usage

```c
#include <ssss.h>

// initialize lib
if (ssss_initialize(SSSS_HEX_MODE_OFF) == -1)
	return -1;

// create shares (secret, threshold, share_number, security_level)
char **shares = ssss_split("your secret", 4, 7, SSSS_DYNAMIC_SECURITY); // security can be [8..1024] multiple of 8

// restore secret from shares (shares, threshold)
char *secret = ssss_combine(shares, 4);

// example with error handling
char *secret = NULL
if ((secret = ssss_combine(shares, 4)) == NULL) {
	printf("%s\n", ssss_get_error_str());
}

// cleanup
ssss_release();
```