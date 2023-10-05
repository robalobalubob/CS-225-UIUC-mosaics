/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    vector<Point<3>> pointies;
    map<Point<3>, TileImage*> maper;
    for (auto it = theTiles.begin(); it != theTiles.end(); ++it) {
        LUVAPixel pixel = it->getAverageColor();
        Point<3> point = convertToXYZ(pixel);
        pointies.push_back(point);
        maper[point] = &*it;
    }
    KDTree<3> avgcolor(pointies);
    MosaicCanvas * canvas = new MosaicCanvas(theSource.getRows(), theSource.getColumns());
    
    for (int i =0; i < canvas->getRows(); i++) {
        for (int j = 0; j < canvas->getColumns(); j++) {
            Point<3> goal = convertToXYZ(theSource.getRegionColor(i,j));
            Point<3> best = avgcolor.findNearestNeighbor(goal);
            TileImage* corTile = maper[best];
            canvas->setTile(i, j, corTile);
        }
    }
    return canvas;
}

