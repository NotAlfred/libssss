/*
 *  sssslib  -  Copyright 2019 NotAlfred
 */

#ifndef SSSS_SSSS_H
#define SSSS_SSSS_H

#include <csprng.h>

#define SSSS_DYNAMIC_SECURITY 0

#define SSSS_HEX_MODE_OFF 0
#define SSSS_HEX_MODE_ON 1

#define SSSS_INVALID_SYNTAX 1
#define SSSS_INPUT_TOO_LONG 2
#define SSSS_INVALID_SECURITY_LEVEL 3
#define SSSS_ILLEGAL_SHARE_LENGTH 4
#define SSSS_SHARE_SECURITY_LEVEL 5
#define SSSS_INVALID_SHARE 6
#define SSSS_INCONSISTENT_SHARE 7

#define VERSION "0.0.1"

#define MAXDEGREE 1024
#define MAXTOKENLEN 128
#define MAXLINELEN (MAXTOKENLEN + 1 + 10 + 1 + MAXDEGREE / 4 + 10)

/* calculate the secret from a set of shares solving a linear equation system */
#define MPZ_SWAP(A, B) \
  do { mpz_set(h, A); mpz_set(A, B); mpz_set(B, h); } while (0)

#define mpz_lshift(A, B, l) mpz_mul_2exp(A, B, l)
#define mpz_sizeinbits(A) (mpz_cmp_ui(A, 0) ? mpz_sizeinbase(A, 2) : 0)

/* coefficients of some irreducible polynomials over GF(2) */
extern const unsigned char ssss_irreductible_coeff[];

/* hexadecimal mode */
extern int ssss_hex_mode;

/* last error code */
extern int ssss_error;

/* random number generator */
extern CSPRNG ssss_rng;

enum encdec {
    ENCODE,
    DECODE
};

int ssss_initialize(int hex_mode);
void ssss_release(void);

char *ssss_get_error_str();

/* return needs to be free'd */
char **ssss_split(char *secret, int threshold, int shares, int level);

/* return need to be free'd */
char *ssss_combine(char **shares, int threshold);

#endif // SSSS_SSSS_H
