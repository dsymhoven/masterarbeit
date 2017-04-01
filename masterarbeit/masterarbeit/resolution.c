//
//  resolution.c
//  masterarbeit
//
//  Created by David Symhoven on 01.04.17.
//  Copyright © 2017 David Symhoven. All rights reserved.
//

#include "resolution.h"

void initResolution(Resolution *Resolution, double dx, double dy, double dz){
    
    Resolution->dx = dx;
    Resolution->dy = dx;
    Resolution->dz = dx;
}
