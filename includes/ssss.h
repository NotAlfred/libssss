/*
 *  sssslib  -  Copyright 2019 NotAlfred
 */

#ifndef SSSS_SSSS_H
# define SSSS_SSSS_H

void split(char **shares, int threshold, int level);
void combine(char **shares, int threshold); // hex - diffusion

#endif // SSSS_SSSS_H
