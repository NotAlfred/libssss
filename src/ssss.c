/*
 *  sssslib  -  Copyright 2019 NotAlfred
 *  ssss version 0.5  -  Copyright 2005,2006 B. Poettering - http://point-at-infinity.org/ssss/
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 *  02111-1307 USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <errno.h>
#include <sys/mman.h>

#include <csprng.h>
#include <gmp.h>

#include "ssss.h"

const unsigned char ssss_irreductible_coeff[] = {
        4,3,1,5,3,1,4,3,1,7,3,2,5,4,3,5,3,2,7,4,2,4,3,1,10,9,3,9,4,2,7,6,2,10,9,
        6,4,3,1,5,4,3,4,3,1,7,2,1,5,3,2,7,4,2,6,3,2,5,3,2,15,3,2,11,3,2,9,8,7,7,
        2,1,5,3,2,9,3,1,7,3,1,9,8,3,9,4,2,8,5,3,15,14,10,10,5,2,9,6,2,9,3,2,9,5,
        2,11,10,1,7,3,2,11,2,1,9,7,4,4,3,1,8,3,1,7,4,1,7,2,1,13,11,6,5,3,2,7,3,2,
        8,7,5,12,3,2,13,10,6,5,3,2,5,3,2,9,5,2,9,7,2,13,4,3,4,3,1,11,6,4,18,9,6,
        19,18,13,11,3,2,15,9,6,4,3,1,16,5,2,15,14,6,8,5,2,15,11,2,11,6,2,7,5,3,8,
        3,1,19,16,9,11,9,6,15,7,6,13,4,3,14,13,3,13,6,3,9,5,2,19,13,6,19,10,3,11,
        6,5,9,2,1,14,3,2,13,3,1,7,5,4,11,9,8,11,6,5,23,16,9,19,14,6,23,10,2,8,3,
        2,5,4,3,9,6,4,4,3,2,13,8,6,13,11,1,13,10,3,11,6,5,19,17,4,15,14,7,13,9,6,
        9,7,3,9,7,1,14,3,2,11,8,2,11,6,4,13,5,2,11,5,1,11,4,1,19,10,3,21,10,6,13,
        3,1,15,7,5,19,18,10,7,5,3,12,7,2,7,5,1,14,9,6,10,3,2,15,13,12,12,11,9,16,
        9,7,12,9,3,9,5,2,17,10,6,24,9,3,17,15,13,5,4,3,19,17,8,15,6,3,19,6,1
};

int ssss_hex_mode = 0;

int ssss_error = 0;

CSPRNG ssss_rng;

unsigned int degree;
mpz_t poly;

void warning(char *msg)
{
    fprintf(stderr, "WARNING: %s.\n", msg);
}

/* field arithmetic routines */
int field_size_valid(int deg)
{
    return (deg >= 8) && (deg <= MAXDEGREE) && (deg % 8 == 0);
}

/* initialize 'poly' to a bitfield representing the coefficients of an
   irreducible polynomial of degree 'deg' */
void field_init(int deg)
{
    assert(field_size_valid(deg));
    mpz_init_set_ui(poly, 0);
    mpz_setbit(poly, deg);
    mpz_setbit(poly, ssss_irreductible_coeff[3 * (deg / 8 - 1) + 0]);
    mpz_setbit(poly, ssss_irreductible_coeff[3 * (deg / 8 - 1) + 1]);
    mpz_setbit(poly, ssss_irreductible_coeff[3 * (deg / 8 - 1) + 2]);
    mpz_setbit(poly, 0);
    degree = deg;
}

void field_deinit(void)
{
    mpz_clear(poly);
}

/* I/O routines for GF(2^deg) field elements */
int field_import(mpz_t x, const char *s, int hex_mode)
{
    if (hex_mode) {
        if (strlen(s) > degree / 4)
            return (ssss_error = SSSS_INPUT_TOO_LONG);
        if (strlen(s) < degree / 4)
            warning("input string too short, adding null padding on the left");
        if (mpz_set_str(x, s, 16) || (mpz_cmp_ui(x, 0) < 0))
            return (ssss_error = SSSS_INVALID_SYNTAX);
    }
    else {
        int warn = 0;
        if (strlen(s) > degree / 8)
            return (ssss_error = SSSS_INPUT_TOO_LONG);
        for (int i = strlen(s) - 1; i >= 0; i--)
            warn = warn || (s[i] < 32) || (s[i] >= 127);
        if (warn)
            warning("binary data detected, use -x mode instead");
        mpz_import(x, strlen(s), 1, 1, 0, 0, s);
    }
    return 0;
}

char *field_to_string(char *prefix, const mpz_t x, int hex_mode)
{
    char *share = NULL;

    if (hex_mode) {
        char padding[degree / 4];
        int j = 0;
        for (int i = degree / 4 - mpz_sizeinbase(x, 16); i; --i, ++j)
            padding[j] = '0';
        padding[j] = '\0';
        gmp_asprintf(&share, "%s%s%Zx", prefix, padding, x);
    }
    else {
        char buf[MAXDEGREE / 8 + 1] = { 0 };
        size_t t;
        int printable, warn = 0;
        mpz_export(buf, &t, 1, 1, 0, 0, x);
        for (unsigned int i = 0; i < t; ++i) {
            printable = (buf[i] >= 32) && (buf[i] < 127);
            warn = warn || !printable;
            buf[i] = printable ? buf[i] : '.';
        }
        gmp_asprintf(&share, "%s%s", prefix, buf);
        if (warn)
            warning("binary data detected, use -x mode instead");
    }
    return share;
}

char *make_share_str(int prefix, const mpz_t field, int hex_mode)
{
    char prefix_str[13] = { 0 };
    sprintf(prefix_str, "%d-", prefix);
    return field_to_string(prefix_str, field, hex_mode);
}

char *make_secret_str(const mpz_t field, int hex_mode)
{
    return field_to_string("", field, hex_mode);
}

/* basic field arithmetic in GF(2^deg) */
void field_add(mpz_t z, const mpz_t x, const mpz_t y)
{
    mpz_xor(z, x, y);
}

void field_mult(mpz_t z, const mpz_t x, const mpz_t y)
{
    mpz_t b;
    assert(z != y);
    mpz_init_set(b, x);
    if (mpz_tstbit(y, 0))
        mpz_set(z, b);
    else
        mpz_set_ui(z, 0);
    for (unsigned int i = 1; i < degree; i++) {
        mpz_lshift(b, b, 1);
        if (mpz_tstbit(b, degree))
            mpz_xor(b, b, poly);
        if (mpz_tstbit(y, i))
            mpz_xor(z, z, b);
    }
    mpz_clear(b);
}

void field_invert(mpz_t z, const mpz_t x)
{
    mpz_t u, v, g, h;
    assert(mpz_cmp_ui(x, 0));
    mpz_init_set(u, x);
    mpz_init_set(v, poly);
    mpz_init_set_ui(g, 0);
    mpz_set_ui(z, 1);
    mpz_init(h);
    while (mpz_cmp_ui(u, 1)) {
        int i = mpz_sizeinbits(u) - mpz_sizeinbits(v);
        if (i < 0) {
            mpz_swap(u, v);
            mpz_swap(z, g);
            i = -i;
        }
        mpz_lshift(h, v, i);
        mpz_xor(u, u, h);
        mpz_lshift(h, g, i);
        mpz_xor(z, z, h);
    }
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(g);
    mpz_clear(h);
}

/* a 64 bit pseudo random permutation (based on the XTEA cipher) */
void encipher_block(uint32_t *v)
{
    uint32_t sum = 0, delta = 0x9E3779B9;
    int i;
    for (i = 0; i < 32; i++) {
        v[0] += (((v[1] << 4) ^ (v[1] >> 5)) + v[1]) ^ sum;
        sum += delta;
        v[1] += (((v[0] << 4) ^ (v[0] >> 5)) + v[0]) ^ sum;
    }
}

void decipher_block(uint32_t *v)
{
    uint32_t sum = 0xC6EF3720, delta = 0x9E3779B9;
    int i;
    for (i = 0; i < 32; i++) {
        v[1] -= ((v[0] << 4 ^ v[0] >> 5) + v[0]) ^ sum;
        sum -= delta;
        v[0] -= ((v[1] << 4 ^ v[1] >> 5) + v[1]) ^ sum;
    }
}

void encode_slice(uint8_t *data, int idx, int len, void (*process_block)(uint32_t *))
{
    uint32_t v[2];
    int i;
    for (i = 0; i < 2; i++)
        v[i] = data[(idx + 4 * i) % len] << 24 |
               data[(idx + 4 * i + 1) % len] << 16 |
               data[(idx + 4 * i + 2) % len] << 8 | data[(idx + 4 * i + 3) % len];
    process_block(v);
    for (i = 0; i < 2; i++) {
        data[(idx + 4 * i + 0) % len] = v[i] >> 24;
        data[(idx + 4 * i + 1) % len] = (v[i] >> 16) & 0xff;
        data[(idx + 4 * i + 2) % len] = (v[i] >> 8) & 0xff;
        data[(idx + 4 * i + 3) % len] = v[i] & 0xff;
    }
}

void encode_mpz(mpz_t x, enum encdec mode)
{
    uint8_t v[(MAXDEGREE + 8) / 16 * 2];
    size_t t;

    memset(v, 0, (degree + 8) / 16 * 2);
    mpz_export(v, &t, -1, 2, 1, 0, x);
    if (degree % 16 == 8)
        v[degree / 8 - 1] = v[degree / 8];
    if (mode == ENCODE)             /* 40 rounds are more than enough!*/
        for (int i = 0; i < 40 * ((int) degree / 8); i += 2)
            encode_slice(v, i, degree / 8, encipher_block);
    else
        for (int i = 40 * (degree / 8) - 2; i >= 0; i -= 2)
            encode_slice(v, i, degree / 8, decipher_block);
    if (degree % 16 == 8) {
        v[degree / 8] = v[degree / 8 - 1];
        v[degree / 8 - 1] = 0;
    }
    mpz_import(x, (degree + 8) / 16, -1, 2, 1, 0, v);
    assert(mpz_sizeinbits(x) <= degree);
}

/* evaluate polynomials efficiently */
void horner(int n, mpz_t y, const mpz_t x, const mpz_t coeff[])
{
    mpz_set(y, x);
    for (int i = n - 1; i; --i) {
        field_add(y, y, coeff[i]);
        field_mult(y, y, x);
    }
    field_add(y, y, coeff[0]);
}

//todo: there must be a way to make this faster as threshold increases (20seconds for t=200 s=1000)
int restore_secret(int n, void *A, mpz_t b[])
{
    mpz_t (*AA)[n] = (mpz_t (*)[n]) A;
    int found;
    mpz_t h;
    mpz_init(h);
    int j;

    for (int i = 0; i < n; i++) {
        if (!mpz_cmp_ui(AA[i][i], 0)) {
            for (j = i + 1, found = 0; j < n; ++j) {
                if (mpz_cmp_ui(AA[i][j], 0)) {
                    found = 1;
                    break;
                }
            }
            if (!found)
                return -1;
            for (int k = i; k < n; ++k)
                MPZ_SWAP(AA[k][i], AA[k][j]);
            MPZ_SWAP(b[i], b[j]);
        }
        for (j = i + 1; j < n; ++j) {
            if (mpz_cmp_ui(AA[i][j], 0)) {
                for (int k = i + 1; k < n; ++k) {
                    field_mult(h, AA[k][i], AA[i][j]);
                    field_mult(AA[k][j], AA[k][j], AA[i][i]);
                    field_add(AA[k][j], AA[k][j], h);
                }
                field_mult(h, b[i], AA[i][j]);
                field_mult(b[j], b[j], AA[i][i]);
                field_add(b[j], b[j], h);
            }
        }
    }
    field_invert(h, AA[n - 1][n - 1]);
    field_mult(b[n - 1], b[n - 1], h);
    mpz_clear(h);
    return 0;
}

/* routines for the random number generator */
void rng_generate(mpz_t x)
{
    char buf[MAXDEGREE / 8];
    csprng_get(ssss_rng, buf, degree / 8); //todo: error handling
    mpz_import(x, degree / 8, 1, 1, 0, 0, buf);
}

int rng_init(void)
{
    ssss_rng = csprng_create(ssss_rng);
    if (!ssss_rng)
        return -1;
    return 0;
}

void rng_release(void)
{
    ssss_rng = csprng_destroy(ssss_rng);
}

void try_lock_mem(void)
{
#if !NOMLOCK
    if (mlockall(MCL_CURRENT | MCL_FUTURE) < 0)
        warning(strerror(errno));
#endif
}

void try_unlock_mem(void)
{
#if !NOMLOCK
    if (munlockall() == -1)
        warning(strerror(errno));
#endif
}

char *ssss_get_error_str()
{
    static char error[256] = { 0 };

    switch (ssss_error) {
        case 0:
            snprintf(error, 18, "No error detected");
            break;
        case SSSS_INVALID_SYNTAX:
            snprintf(error, 15, "Invalid syntax");
            break;
        case SSSS_INPUT_TOO_LONG:
            snprintf(error, 22, "Input string too long");
            break;
        case SSSS_INVALID_SECURITY_LEVEL:
            snprintf(error, 43, "Security level invalid (secret too long?)");
            break;
        case SSSS_ILLEGAL_SHARE_LENGTH:
            snprintf(error, 26, "Share has illegal length");
            break;
        case SSSS_SHARE_SECURITY_LEVEL:
            snprintf(error, 39, "Shares have different security levels");
            break;
        case SSSS_INVALID_SHARE:
            snprintf(error, 14, "Invalid share");
            break;
        case SSSS_INCONSISTENT_SHARE:
            snprintf(error, 60, "Shares inconsistent. Perhaps a single share was used twice");
            break;
    }
    ssss_error = 0;
    return error;
}

/* Prompt for a secret, generate shares for it */
//todo: error mgmt on params
char **ssss_split(char *secret, int threshold, int share_number, int security_level)
{
    char **shares;
    //unsigned int fmt_len;
    mpz_t x, y, coeff[threshold];
    for (int i = share_number, fmt_len = 1; i >= 10; i /= 10, ++fmt_len);

    if (strlen(secret) >= MAXLINELEN)
        return NULL;

    if (security_level == SSSS_DYNAMIC_SECURITY) {
        security_level = ssss_hex_mode ? 4 * ((strlen(secret) + 1) & ~1) : 8 * strlen(secret);
        if (!field_size_valid(security_level)) {
            ssss_error = SSSS_INVALID_SECURITY_LEVEL;
            return NULL;
        }
    }

    field_init(security_level);

    mpz_init(coeff[0]);
    if (field_import(coeff[0], secret, ssss_hex_mode) != 0)
        return NULL;

    // diffusion layer
    if (degree >= 64)
        encode_mpz(coeff[0], ENCODE);
    else
        warning("security level too small for the diffusion layer");

    for (int i = 1; i < threshold; ++i) {
        mpz_init(coeff[i]);
        rng_generate(coeff[i]);
    }

    shares = calloc(share_number, sizeof(char *));

    mpz_init(x);
    mpz_init(y);
    for (int i = 0; i < share_number; ++i) {
        mpz_set_ui(x, i + 1);
        horner(threshold, y, x, (const mpz_t *) coeff);
        shares[i] = make_share_str(i + 1, y, 1);
        //if (opt_token) //todo: possibility to add a prefix token
          //  printf("%s-", opt_token);
    }
    mpz_clear(x);
    mpz_clear(y);

    for (int i = 0; i < threshold; ++i)
        mpz_clear(coeff[i]);
    field_deinit();
    return shares;
}

/* Prompt for shares, calculate the secret */
//todo: error mgmt on params
char *ssss_combine(char **shares, int threshold)
{
    mpz_t A[threshold][threshold], y[threshold], x;
    char *a, *b;
    unsigned s = 0;

    mpz_init(x);
    for (int i = 0; i < threshold; ++i) {
        if (!(a = strchr(shares[i], '-'))) {
            ssss_error = SSSS_INVALID_SYNTAX;
            return NULL;
        }
        *a++ = 0;
        if ((b = strchr(a, '-')))
            *b++ = 0;
        else
            b = a, a = shares[i];

        if (!s) {
            s = 4 * strlen(b);
            if (!field_size_valid(s)) {
                ssss_error = SSSS_ILLEGAL_SHARE_LENGTH;
                return NULL;
            }
            field_init(s);
        }

        else if (s != 4 * strlen(b)) {
            ssss_error = SSSS_SHARE_SECURITY_LEVEL;
            return NULL;
        }

        int prefix;
        if (!(prefix = atoi(a))) {
            ssss_error = SSSS_INVALID_SHARE;
            return NULL;
        }
        mpz_set_ui(x, prefix);
        mpz_init_set_ui(A[threshold - 1][i], 1);
        for (int j = threshold - 2; j >= 0; --j) {
            mpz_init(A[j][i]);
            field_mult(A[j][i], A[j + 1][i], x);
        }
        mpz_init(y[i]);
        if (field_import(y[i], b, 1) != 0)
            return NULL;
        field_mult(x, x, A[0][i]);
        field_add(y[i], y[i], x);
    }
    mpz_clear(x);

    if (restore_secret(threshold, A, y)) {
        ssss_error = SSSS_INCONSISTENT_SHARE;
        return NULL;
    }

    // diffusion layer
    if (degree >= 64)
        encode_mpz(y[threshold - 1], DECODE);
    else
        warning("security level too small for the diffusion layer");

    char *share = make_secret_str(y[threshold - 1], ssss_hex_mode);

    for (int i = 0; i < threshold; ++i) {
        for (int j = 0; j < threshold; ++j)
            mpz_clear(A[i][j]);
        mpz_clear(y[i]);
    }
    field_deinit();
    return share;
}

int ssss_initialize(int hex_mode)
{
    ssss_hex_mode = hex_mode;
    ssss_error = 0;
    if (rng_init() == -1)
        return -1;
    try_lock_mem();
    return 0;
}

void ssss_release(void)
{
    rng_release();
    try_unlock_mem();
}

int main(int argc, char *argv[])
{
    if (ssss_initialize(SSSS_HEX_MODE_OFF) == -1)
        return -1;
    char **shares = ssss_split("totototototototafzefzv", 200, 1000, SSSS_DYNAMIC_SECURITY); // [8..1024] multiple of 8
    for (int i = 0; i < 1000; ++i) {
        printf("%s\n", shares[i]);
    }
    char *share = ssss_combine(shares, 200);
    printf("%s\n", share);
    ssss_release();
    return 0;
}
