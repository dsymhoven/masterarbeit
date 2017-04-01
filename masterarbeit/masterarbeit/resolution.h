//
//  resolution.h
//  masterarbeit
//
//  Created by David Symhoven on 01.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#ifndef resolution_h
#define resolution_h

#include <stdio.h>

#endif /* resolution_h */

struct Resolution{
    double dx;
    double dy;
    double dz;
};

typedef struct Resolution Resolution;

void initResolution(Resolution *Resolution, double dx, double dy, double dz);
