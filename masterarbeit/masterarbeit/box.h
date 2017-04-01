//
//  box.h
//  masterarbeit
//
//  Created by David Symhoven on 01.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#ifndef box_h
#define box_h

#include <stdio.h>

#endif /* box_h */

struct Box{
    int numberOfGridPointsInX;
    int numberOfGridPointsInY;
    int numberOfGridPointsInZ;
};

typedef struct Box Box;

void initBox(Box *Box, const int numberOfGridPointsInX, const int numberOfGridPointsInY, const int numberOfGridPointsInZ);
