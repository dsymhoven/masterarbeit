//
//  experiment.c
//  masterarbeit
//
//  Created by David Symhoven on 02.11.16.
//  Copyright © 2016 David Symhoven. All rights reserved.
//

#include "experiment.h"
#include "stdlib.h"
#include "stdio.h"
#include "grid.h"


void experiment1(){
    Grid Grid;
    initGrid(&Grid, 256, 32);
    
    
    printf("%d\n", Grid.numberOfGridPoints);
    freeMemoryOn(&Grid);
}
