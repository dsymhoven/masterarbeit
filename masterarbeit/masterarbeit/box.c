//
//  box.c
//  masterarbeit
//
//  Created by David Symhoven on 01.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#include "box.h"

void initBox(Box *Box, const int numberOfGridPointsInX, const int numberOfGridPointsInY, const int numberOfGridPointsInZ){
    Box->numberOfGridPointsInX = numberOfGridPointsInX;
    Box->numberOfGridPointsInY = numberOfGridPointsInY;
    Box->numberOfGridPointsInZ = numberOfGridPointsInZ;
}
