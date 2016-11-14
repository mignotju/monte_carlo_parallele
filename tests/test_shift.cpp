/* 
 * File:   test_shift.cpp
 * Author: gononl
 *
 * Created on September 21, 2016, 2:05 PM
 */

#include <cstdlib>
#include <iostream>

#include "../src/BlackScholesModel.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    PnlMat *path = pnl_mat_create_from_scalar(10,6,10);
    pnl_mat_print(path);
    cout << " " << endl;
    BlackScholesModel *bs = new BlackScholesModel();
    
    PnlMat *shiftPath = pnl_mat_create(10,6);
    
    bs->shiftAsset(shiftPath, path, 3, -1, 10, 1);

    pnl_mat_print(shiftPath);
    
    return 0;
}

