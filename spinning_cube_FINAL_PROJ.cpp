//Justin Brown FINAL PROJECT c++
// 9/24/2022 ->
//Display a spinning geometric cube in the terminal
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cstring>
#include <string>
#include <chrono>
#include <thread>
#include <array>
#define print_debug_information false

//Structs
struct geoSolid{
  double* points;
  int numPoints;
  int* projection;
  int* connectedPoints;
  int numConnectedPoints;
  int* faces;
  int numFaces;
  int* facesToPrint;
  int numFacesToPrint;
  double* lines;
};
struct perspectivePoint {
    int planeIndex, coord1, coord2;
    double* points;
    double planeVal;
};
struct displayConfig {
    double xyRatio = 0.5;
    int gridWidth, gridHeight;
    int max_screen_width = 300;
    int max_screen_height = 200;
    double speedScalar = 7;
    int delayNanos = 10000 + 7500 / speedScalar;
    double speedPiDivisor = (1000 / speedScalar) * 2;
    double thetaX = M_PI / speedPiDivisor;
    double thetaY = M_PI / speedPiDivisor;
    double thetaZ = M_PI / speedPiDivisor;
    double pointScalar = 40;
    double zoomScalar = 2;
    double pushUp = 4;
    double pushRight = 0;
    double xPower = 0.2;
    double zPower = 0.5;
    double yPower = 1.1;
    char background_char = ' ';
    char connectionChars[3] = { '+','%','#' };
    char pointChars[3] = { '*', 'o', '@' };
};
// =================prototypes======================
//most important functions
void rotate(geoSolid& shape, double thetaX, double thetaY, double thetaZ);
void project(geoSolid& shape, perspectivePoint(&perspectives)[2], displayConfig& dispConf);
void getCorrectFaces(geoSolid& shape, std::string* canvas, displayConfig& dispConf);
int* paint(std::string* canvas, geoSolid& shape, displayConfig& dispConf, bool doShadow);
//point realignment functions
void functionalize(geoSolid& shape);
void fixSideLength(geoSolid& shape, double targetLineLength);
void vectorSpread(geoSolid& shape, int p1Index, int p2Index, double targetLineLength);
void fixRadialDistance(geoSolid& shape, double targetDistance);
void scale_points(geoSolid& shape, double scaler);
//helper functions
void arr_fill(double* target, double* src, int rows, int cols);
void arr_fill(int* target, int* src, int rows, int cols);
double distance(double* p1, double* p2);
void addFaceToPrint(geoSolid& shape, int face);
bool inBounds(double num, double bound1, double bound2);
int* indexOfIntArr(int* elem, int* arr, int rowlen, int collen, int elemlen);
int indexOfCharArr(char elem, char* arr, int len);
bool checkConnectionGoodToPrint(geoSolid& shape, int connectionIndex);
int* get_screen_dimensions(int maxWidth, int maxHeight);
//display functions
void print_startup_banner();
void print(std::string* canvas, displayConfig& dispConf, int lowestY);
void fill_newlines(int count);
void print_points(double* pts);
//struct initializers
void init_shape(int shapeNum, geoSolid& shape);
void init_perspectives(perspectivePoint* perspectives);


int main(int argc, char const* argv[]) {
    /*==============SETUP==============*/
    double runningthetaX;
    double runningthetaY;
    double runningthetaZ;
    double center[3] = { 0,0,0 };
    double targetLineLength = 0;
    double targetRadialDistance = 0;
    int loop_max_iterations = 1000000;
    int* dimensions = nullptr;
    displayConfig dispConf;
    perspectivePoint perspectives[2];
    geoSolid shape;
    std::string* canvas;
    
    init_shape(0, shape);
    init_perspectives(perspectives);

    scale_points(shape, dispConf.pointScalar);//Scale points ( 1,1,1 -> 5,5,5 ... etc)
    targetLineLength = distance(shape.points + (*(shape.connectedPoints) * 3), (shape.points + (*(shape.connectedPoints + 1) * 3)));
    targetRadialDistance = distance(shape.points, center);

    print_startup_banner(); //Ascii Graphic
    dimensions = get_screen_dimensions(dispConf.max_screen_width, dispConf.max_screen_height);//Get Dimensions From User & set relevent vars
    dispConf.gridWidth = dimensions[0];
    dispConf.gridHeight = dimensions[1];
    canvas = new std::string[dispConf.gridHeight];
    int lowestY = 0;
    int highestY = dispConf.gridHeight;

    fill_newlines(500); //Fills Terminal with newline chars to avoid jittering effect

    if (print_debug_information) printf("X:%i  Y:%i\n", dispConf.gridWidth, dispConf.gridHeight); 
    /*==============MAIN LOOP==============*/
    auto startTime = std::chrono::system_clock::now();
    double ellapsed;
    for (int i = 0; i < loop_max_iterations; i++) {
        //Calculate thetas to use this iteration
        ellapsed = (std::chrono::system_clock::now() - startTime).count();
        runningthetaX = dispConf.thetaX * (std::cos(ellapsed * dispConf.xPower) + 13 / 12) * dispConf.xPower;
        runningthetaY = dispConf.thetaY * (std::sin(ellapsed * dispConf.yPower) + 1) * dispConf.yPower;
        runningthetaZ = dispConf.thetaZ * (std::sin(ellapsed * dispConf.zPower) + 1) * dispConf.zPower;
        if (i > 50) {
            if (print_debug_information) printf("\nThetas- X:%f Y:%f Z:%f |", runningthetaX, runningthetaY, runningthetaZ);
            if (print_debug_information) printf(" YRange: %i -> %i |", lowestY, highestY);
        }
        //rotate points
        rotate(shape, runningthetaX, runningthetaY, runningthetaZ);

        //fix misalignments in points due to floating point arithmetic 
        fixSideLength(shape, targetLineLength);
        fixRadialDistance(shape, targetRadialDistance);

        //project 3d points to 2d planes 
        project(shape, perspectives, dispConf);

        //calculate y=mx+b type equations for each projected connection
        functionalize(shape);

        //use intersections of connections and z coordinates of faces to determine correct faces to display
        getCorrectFaces(shape, canvas, dispConf);

        //paint lines on canvas, lowestY/HighestY likely unnessesary, currently unused
        int* yRangePointer = paint(canvas, shape, dispConf, false);
        lowestY = *(yRangePointer);
        highestY = *(yRangePointer + 1);
        //sleep
        if (i > 50) { // let the loop run for a bit so the yRange stablizes
            std::this_thread::sleep_for(std::chrono::nanoseconds(dispConf.delayNanos));

            //print canvas
            print(canvas, dispConf, lowestY);
        }
        


    }
    return 0;
}


/*######################## Main Functions ########################*/
void rotate(geoSolid& shape, double thetaX, double thetaY, double thetaZ) {
    //Rotate points about the origin according to each directions theta value

    //Calculate sin & cos for each theta value
    //X
    double sinThetaX = std::sin(thetaX);
    double cosThetaX = std::cos(thetaX);
    //Y
    double sinThetaY = std::sin(thetaY);
    double cosThetaY = std::cos(thetaY);
    //Z
    double sinThetaZ = std::sin(thetaZ);
    double cosThetaZ = std::cos(thetaZ);

    for (int i = 0;i < shape.numPoints;i++) //rotate each point
    {
        double x = *(shape.points + (i * 3) + 0);
        double y = *(shape.points + (i * 3) + 1);
        double z = *(shape.points + (i * 3) + 2);
        //X
        y = (y * cosThetaX - z * sinThetaX);
        z = (z * cosThetaX + y * sinThetaX);
        //Y
        x = (x * cosThetaY - z * sinThetaY);
        z = (z * cosThetaY + x * sinThetaY);
        //Z
        x = (x * cosThetaZ - y * sinThetaZ);
        y = (y * cosThetaZ + x * sinThetaZ);

        //apply rotation to class
        *(shape.points + (i * 3) + 0) = x;
        *(shape.points + (i * 3) + 1) = y;
        *(shape.points + (i * 3) + 2) = z;
    }
}

void project(geoSolid& shape, perspectivePoint(&perspectives)[2], displayConfig& dispConf) {
    /* Project 3d points objPoints onto the 2d plane on z=planeZ, from the perspective of the point eyePoints  */
    perspectivePoint sun = perspectives[0]; //for shadow
    perspectivePoint eye = perspectives[1]; //where the viewer is relative to the shape
    int* windowPoints = shape.projection + (2 * shape.numPoints * 3); //point -> viewer
    int* shadowPoints = shape.projection + (1 * shape.numPoints * 3); //point -> shadow
    int* shadowWindowPoints = shape.projection; //shadow -> viewer

    for (int i = 0; i < shape.numPoints; i++) {
        // #### point -> viewer #### ================================================================================================================
        double x = *(shape.points + (i * 3) + 0);
        double y = *(shape.points + (i * 3) + 1);
        double z = *(shape.points + (i * 3) + 2);
        double coords[3] = { x,y,z };
        //calculation
        double projx = (((eye.points[0] - coords[0]) / (eye.points[2] - coords[2])) * eye.planeVal + eye.points[0]) * dispConf.zoomScalar + dispConf.pushRight;
        double projy = (((eye.points[1] - coords[1]) / (eye.points[2] - coords[2])) * eye.planeVal + eye.points[1]) * dispConf.zoomScalar * dispConf.xyRatio;
        //adding 1/2 so int cast rounds to nearest integer
        projx += 0.5;
        projy += 0.5;
        //apply calculation to projection array
        *(windowPoints + (i * 3) + 0) = static_cast<int>(projx);
        *(windowPoints + (i * 3) + 1) = static_cast<int>(projy);
        *(windowPoints + (i * 3) + 2) = static_cast<int>(eye.planeIndex + 0.5);

        // #### point -> shadow #### ================================================================================================================
        //calculation
        projx = (((sun.points[0] - coords[0]) / (sun.points[1] - coords[1])) * sun.planeVal + sun.points[0]);
        double projz = (((sun.points[2] - coords[2]) / (sun.points[1] - coords[1])) * sun.planeVal + sun.points[2]);
        //adding 1/2 so int cast rounds to nearest integer
        projx += 0.5;
        projz += 0.5;
        //apply calculation to projection array
        *(shadowPoints + (i * 3) + 0) = static_cast<int>(projx);
        *(shadowPoints + (i * 3) + 2) = static_cast<int>(projz);
        *(shadowPoints + (i * 3) + 1) = static_cast<int>(sun.planeVal);

        // #### shadow -> viewer #### ================================================================================================================
        coords[0] = *(shadowPoints + (i * 3) + 0);
        coords[1] = *(shadowPoints + (i * 3) + 1);
        coords[2] = *(shadowPoints + (i * 3) + 2);
        //calculation
        projx = (((eye.points[0] - coords[0]) / (eye.points[2] - coords[2])) * eye.planeVal + eye.points[0]) * dispConf.zoomScalar + dispConf.pushRight;
        projy = (((eye.points[1] - coords[1]) / (eye.points[2] - coords[2])) * eye.planeVal + eye.points[1]) * dispConf.zoomScalar * dispConf.xyRatio;
        //adding 1/2 so int cast rounds to nearest integer
        projx += 0.5;
        projy += 0.5;
        //apply calculation to projection array
        *(shadowWindowPoints + (i * 3) + 0) = static_cast<int>(projx);
        *(shadowWindowPoints + (i * 3) + 1) = static_cast<int>(projy);
        *(shadowWindowPoints + (i * 3) + 2) = static_cast<int>(coords[2] + 0.5);
    }
}

void getCorrectFaces(geoSolid& shape, std::string* canvas, displayConfig& dispConf) {
    //clears canvas because this function is the first place where it might be written to
    static int* lastFacesToPrint = new int[6];
    for (int i = 0; i < dispConf.gridHeight; i++) {
        canvas[i] = std::string(dispConf.gridWidth, dispConf.background_char);
    }
    //Looks for connection intersections by using the connections in y= mx +b | x > x1 and x < x2 calculated by the functionalize function
    int yAxis = (static_cast<double>(dispConf.gridWidth) + 0.5) / 2;
    int xAxis = (static_cast<double>(dispConf.gridHeight) + 0.5) / 2;
    shape.numFacesToPrint = 0;
    int numFound = 0;
    int foundIntersection[12] = {};
    for (int i = 0; i < shape.numConnectedPoints; i++) {
        //connection 1, looking for another connection that intersects it
        double m1 = *(shape.lines + (i * 4));
        double b1 = *(shape.lines + (i * 4) + 1);
        double startX1 = *(shape.lines + (i * 4) + 2);
        double endX1 = *(shape.lines + (i * 4) + 3);
        double startY1 = m1 * startX1 + b1;
        double endY1 = m1 * endX1 + b1;
        for (int w = 0; w < shape.numConnectedPoints; w++) {
            //connection 2, check if it intersects
            double m2 = *(shape.lines + (w * 4));
            double b2 = *(shape.lines + (w * 4) + 1);
            double startX2 = *(shape.lines + (w * 4) + 2);
            double endX2 = *(shape.lines + (w * 4) + 3);
            bool isNotSame = w != i; //make sure its not the same line
            if (isNotSame) {
                double startY2 = m2 * startX2 + b2;
                double endY2 = m2 * endX2 + b2;
                double x = (b1 - b2) / (m2 - m1);
                double y;
                y = x * m1 + b1;
                //Do bounds checking to avoid subscript error
                if (x + yAxis >= 0 && x + yAxis < yAxis * 2 && xAxis - y + 0.5 + dispConf.pushUp + 7 >= 0 && xAxis - y + 0.5 + dispConf.pushUp + 7 < xAxis * 2) {
                    //Do bounds checking to make sure intersection is within the lines, not outside the shape
                    if (inBounds(x, startX1, endX1) && inBounds(x, startX2, endX2))
                    {
                        // intersection coordinates found, get indecies of the connections
                        int c1i1 = *(shape.connectedPoints + (i * 2));
                        int c1i2 = *(shape.connectedPoints + (i * 2) + 1);
                        int c2i1 = *(shape.connectedPoints + (w * 2));
                        int c2i2 = *(shape.connectedPoints + (w * 2) + 1);

                        //ensure none of the connections share a point, this avoids detecting corners as intersections
                        if (c1i1 != c2i2 && c1i1 != c2i1 && c1i2 != c2i1 && c1i2 != c2i2) {
                            //record this so we dont repeat the intersection
                            foundIntersection[i] = w;
                            //calculate the average z coords of each connection to find which one is 'on top' i.e. should be displayed 
                            double c1avgZ = (*(shape.points + (3 * c1i1) + 2) + *(shape.points + (3 * c1i2) + 2)) / 2;
                            double c2avgZ = (*(shape.points + (3 * c2i1) + 2) + *(shape.points + (3 * c2i2) + 2)) / 2;
                            
                            int* facesToPaint = new int[2];
                            if (c1avgZ - c2avgZ > 8) { //if connection1 is on top
                                //get faces that have this connection as a line segment
                                facesToPaint = indexOfIntArr((shape.connectedPoints + (i * 2)), shape.faces, 6, 4, 2);
                                
                                //checks to make sure we're not throwing a -1 into the array
                                if (facesToPaint[0] >= 0) {
                                    addFaceToPrint(shape, facesToPaint[0]);
                                }
                                if (facesToPaint[1] >= 0) {
                                    addFaceToPrint(shape, facesToPaint[1]);
                                }
                                //increment numFound
                                numFound++;

                                //if print_debug_information is Enabled, mark the intersection with an X
                                if (print_debug_information) {
                                    printf("Print Faces %i & %i ", facesToPaint[0], facesToPaint[1]);
                                    canvas[static_cast<int>(xAxis - y + 0.5 + dispConf.pushUp)][static_cast<int>(x + 0.5 + yAxis)] = 'X';
                                }
                                
                            }
                            else if (c2avgZ - c1avgZ > 8) {//if connection2 is on top
                                facesToPaint = indexOfIntArr((shape.connectedPoints + (w * 2)), shape.faces, 6, 4, 2);
                                if (facesToPaint[0] >= 0) {
                                    addFaceToPrint(shape, facesToPaint[0]);
                                }
                                if (facesToPaint[1] >= 0) {
                                    addFaceToPrint(shape, facesToPaint[1]);
                                }
                                numFound++;
                                if (print_debug_information) {
                                    printf("Print Faces %i & %i ", facesToPaint[0], facesToPaint[1]);
                                    canvas[static_cast<int>(xAxis - y + 0.5 + dispConf.pushUp)][static_cast<int>(x + 0.5 + yAxis)] = 'X';
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    if (shape.numFacesToPrint == 0) { //Occurs when one face is 'facing' the viewer, no connections intersect on the projection
        double maxZ = -100;
        int maxFaceIndex = -1;
        for (int i = 0; i < shape.numFaces; i++) {//cycle through all the faces & their points z-Values take the face with the highest average z value and display that
            double faceZTotal = 0;
            int* face = shape.faces + (4 * i);
            for (int f = 0; f < 4; f++) {
                double* point = shape.points + (3 * *(face + f));
                faceZTotal += point[2];
            }
            if (faceZTotal > maxZ) {
                maxFaceIndex = i;
                maxZ = faceZTotal;
            }
        }
        shape.numFacesToPrint = 1;
        shape.facesToPrint[0] = maxFaceIndex;
        if (print_debug_information) printf("Print Face %i ", maxFaceIndex);
    }
}

int* paint(std::string* canvas, geoSolid& shape, displayConfig& dispConf, bool doShadow)
{   /*  Convert projected points to printable canvas of strings  */
    int planeIndex = 2;
    if (doShadow) {
        planeIndex = 0;
    }
    //setup variables
    static int highestY = -500;
    static int lowestY = 500;
    int yAxis = (static_cast<double>(dispConf.gridWidth) + 0.5) / 2;
    int xAxis = (static_cast<double>(dispConf.gridHeight) + 0.5) / 2;

    for (int i = 0; i < shape.numConnectedPoints; i++) {//iterate through connected point pairs
        bool test = false;
        if(doShadow || checkConnectionGoodToPrint(shape, i) || test)
        {
            //Two points connected by an edge
            int* p1;
            int* p2;

            if (*(shape.projection + (planeIndex * shape.numPoints * 3) + (*(shape.connectedPoints + (i * 2) + 0) * 3) + 0) > *(shape.projection + (planeIndex * shape.numPoints * 3) + (*(shape.connectedPoints + (i * 2) + 1) * 3) + 0)) {
                p1 = (shape.projection + (planeIndex * shape.numPoints * 3) + (*(shape.connectedPoints + (i * 2) + 1) * 3));
                p2 = (shape.projection + (planeIndex * shape.numPoints * 3) + (*(shape.connectedPoints + (i * 2) + 0) * 3));
            }
            else {
                p1 = (shape.projection + (planeIndex * shape.numPoints * 3) + (*(shape.connectedPoints + (i * 2) + 0) * 3));
                p2 = (shape.projection + (planeIndex * shape.numPoints * 3) + (*(shape.connectedPoints + (i * 2) + 1) * 3));
            }


            //Find Slope
            double xDiff = p1[0] - p2[0];
            double yDiff = p1[1] - p2[1];
            bool vertical = false;
            double slope;
            if (std::abs(xDiff) <= 1) {
                vertical = true;
            }
            else {
                slope = yDiff / xDiff;
                if (std::abs(slope) >= 7) {
                    vertical = true;
                }
                //else {
                //    
                //}
            }
            //Plot points on cavas
            if (vertical) {
                int x = p1[0];
                int x2 = p2[0];
                int y1 = p1[1];
                int y2 = p2[1];
                if (p1[1] > p2[1]) {
                    y2 = p1[1];
                    y1 = p2[1];
                }
                int connectionCharIndex;
                if (p1[2] > 0 && p2[2] > 0) {
                    connectionCharIndex = 2;
                }
                else if (p1[2] > 0 || p2[2] > 0) {
                    connectionCharIndex = 1;
                }
                else {
                    connectionCharIndex = 0;
                }
                //if (y1 >= y2) {
                //    std::cout << "";
                //}
                //if (y1 < 0) {
                //    y1 = 0;
                //}
                //if (y2 > dispConf.gridHeight) {
                //    y2 = dispConf.gridHeight;
                //}
                for (int y = y1; y <= y2; y++)
                {
                    int yOnCanvas = xAxis - y + dispConf.pushUp;
                    if (x + yAxis >= 0 && x + yAxis < dispConf.gridWidth && yOnCanvas >= 0 && yOnCanvas < dispConf.gridHeight) {
                        char thisChar = canvas[yOnCanvas][x + yAxis];
                        bool isNotOverLapping = thisChar == dispConf.background_char || indexOfCharArr(thisChar, dispConf.connectionChars, 3) < connectionCharIndex;
                        bool isNotX = thisChar != 'X';
                        if (isNotOverLapping && isNotX) {
                            canvas[yOnCanvas][x + yAxis] = dispConf.connectionChars[connectionCharIndex];
                            if (yOnCanvas > highestY) {
                                highestY = yOnCanvas;
                            }
                            if (yOnCanvas < lowestY) {
                                lowestY = yOnCanvas;
                            }
                        }
                    }
                }
            }
            else {
                int x1, x2, y1, y2;
                if (p1[0] > p2[0]) {
                    x2 = p1[0];
                    x1 = p2[0];
                    y2 = p1[1];
                    y1 = p2[1];
                }
                else {
                    x1 = p1[0];
                    x2 = p2[0];
                    y1 = p1[1];
                    y2 = p2[1];
                }
                if (x1 + yAxis < 0) {
                    x1 = 0;
                }
                if (x2 + yAxis >= dispConf.gridWidth) {
                    x2 = dispConf.gridWidth - yAxis;
                }
                for (int x = x1; x < x2; x++) {
                    //Choose which char to draw lines with
                    int connectionCharIndex;
                    if (doShadow)
                    {
                        connectionCharIndex = 0;
                    }
                    else if (p1[2] > 0 && p2[2] > 0) {
                        connectionCharIndex = 2;
                    }
                    else if (p1[2] > 0 || p2[2] > 0) {
                        connectionCharIndex = 1;
                    }
                    else {
                        connectionCharIndex = 0;
                    }
                    //int y = static_cast<double>(x - x1) * slope + 0.5;
                    int yStart = xAxis - static_cast<double>(x - x1) * slope + 0.5 - y1 + dispConf.pushUp;
                    int yEnd = xAxis - static_cast<double>(x - x1 + 1) * slope + 0.5 - y1 + dispConf.pushUp;

                    int incDirection = (yStart <= yEnd) ? 1 : -1;

                    for (int yToPlot = yStart; yToPlot != yEnd;yToPlot += incDirection) {
                        if (yToPlot < dispConf.gridHeight && yToPlot >= 0 && x + yAxis < dispConf.gridWidth && x + yAxis >= 0) {
                            char thisChar = canvas[yToPlot][x + yAxis];
                            bool isNotOverLapping = thisChar == dispConf.background_char || indexOfCharArr(thisChar, dispConf.connectionChars, 3) < connectionCharIndex;
                            bool isNotX = thisChar != 'X';
                            if (isNotOverLapping && isNotX) {
                                canvas[yToPlot][x + yAxis] = dispConf.connectionChars[connectionCharIndex];
                                if (yToPlot > highestY) {
                                    highestY = yToPlot;
                                }
                                if (yToPlot < lowestY) {
                                    lowestY = yToPlot;
                                }
                            }
                        }
                        //else {
                        //    printf("Fail Shape X:%i Y:%i", x + yAxis, yToPlot);
                        //}
                    }
                    for (int yToPlot = yStart; yToPlot != yEnd + 1;yToPlot += incDirection) {
                        if (yToPlot < dispConf.gridHeight && yToPlot >= 0 && x + yAxis < dispConf.gridWidth && x + yAxis >= 0) {
                            char thisChar = canvas[yToPlot][x + yAxis];
                            bool isNotOverLapping = thisChar == dispConf.background_char || indexOfCharArr(thisChar, dispConf.connectionChars, 3) < connectionCharIndex;
                            bool isNotX = thisChar != 'X';
                            if (isNotOverLapping && isNotX) {
                                canvas[yToPlot][x + yAxis] = dispConf.connectionChars[connectionCharIndex];
                                if (yToPlot > highestY) {
                                    highestY = yToPlot;
                                }
                                if (yToPlot < lowestY) {
                                    lowestY = yToPlot;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    int *returnData = nullptr;
    returnData = new int[2];
    returnData[0] = lowestY;
    returnData[1] = highestY;

    if (!doShadow) {
        paint(canvas, shape, dispConf, true);
    }
    return returnData;
}


/*######################## Point Realignment Functions ########################*/
void functionalize(geoSolid& shape) {
    //get the y=mx + b representation of each projection connection
    int* cubeProjection = shape.projection + (shape.numPoints * 3 * 2);
    for (int i = 0; i < shape.numConnectedPoints; i++) { //for each connection on the projection of the cube
        int p1Index, p2Index;
        p1Index = *(shape.connectedPoints + (i * 2));
        p2Index = *(shape.connectedPoints + (i * 2) + 1);
        int* p1, * p2;
        p1 = (cubeProjection + (p1Index * 3)); 
        p2 = (cubeProjection + (p2Index * 3));
        //if (p2 < p1) {
        //    p1 = (cubeProjection + (p2Index * 3));
        //    p2 = (cubeProjection + (p1Index * 3));
        //}
        double xDiff = *(p2)-*(p1); //the 'run' of the slope formula
        if (xDiff == 0) { //avoid division by 0
            xDiff = 0.0000000001;
        }
        double yDiff = *(p2 + 1) - *(p1 + 1); //the 'rise' of the slope formula
        double m = yDiff / xDiff; //slope aka 'm'
        double b = *(p1 + 1) - (m * *p1); //y intercept aka 'b'
        *(shape.lines + (i * 4) + 0) = m;
        *(shape.lines + (i * 4) + 1) = b;
        *(shape.lines + (i * 4) + 2) = *p1;
        *(shape.lines + (i * 4) + 3) = *p2;

    }
}

void fixSideLength(geoSolid& shape, double targetLineLength) {
    /*calls vectorSpread on each connected point*/
    for (int i = 0; i < shape.numConnectedPoints; i++) {
        vectorSpread(shape, *(shape.connectedPoints + (i * 2) + 0), *(shape.connectedPoints + (i * 2) + 1), targetLineLength);
    }
    if (shape.numPoints == 8 && print_debug_information) {/*measures line distance for cubes*/
        int oneDifferenceConnections[4][2];
        double oneDiffAvg = 0;
        int oneDiffCount = 0;

        int twoDifferenceConnections[4][2];
        double twoDiffAvg = 0;
        int twoDiffCount = 0;

        int fourDifferenceConnections[4][2];
        double fourDiffAvg = 0;
        int fourDiffCount = 0;

        //the connections & points are ordered in such a way that the difference between the connection indecies is the same for parallel connections, allowing us to take the average of
        //the lengths of the 3 groups of parallel connections
        for (int i = 0; i < shape.numConnectedPoints; i++) {
            double lineDistance = distance(shape.points + (3 * *(shape.connectedPoints + (i * 2) + 1)), shape.points + 3 * *(shape.connectedPoints + (i * 2) + 0));
            int connectionFingerprint = *(shape.connectedPoints + (i * 2) + 1) - *(shape.connectedPoints + (i * 2) + 0);
            if (connectionFingerprint == 1) {
                oneDifferenceConnections[oneDiffCount][0] = *(shape.connectedPoints + (i * 2) + 0);
                oneDifferenceConnections[oneDiffCount][1] = *(shape.connectedPoints + (i * 2) + 1);
                oneDiffAvg = (lineDistance + (oneDiffCount * oneDiffAvg)) / (oneDiffCount + 1);
                oneDiffCount++;
            }
            else if (connectionFingerprint == 2) {
                twoDifferenceConnections[twoDiffCount][0] = *(shape.connectedPoints + (i * 2) + 0);
                twoDifferenceConnections[twoDiffCount][1] = *(shape.connectedPoints + (i * 2) + 1);
                twoDiffAvg = (lineDistance + (twoDiffCount * twoDiffAvg)) / (twoDiffCount + 1);
                twoDiffCount++;
            }
            else {
                fourDifferenceConnections[fourDiffCount][0] = *(shape.connectedPoints + (i * 2) + 0);
                fourDifferenceConnections[fourDiffCount][1] = *(shape.connectedPoints + (i * 2) + 1);
                fourDiffAvg = (lineDistance + (fourDiffCount * fourDiffAvg)) / (fourDiffCount + 1);
                fourDiffCount++;
            }
        }
        printf(" Avg Line Length- X:%f Y:%f Z:%f Target:%f | ", oneDiffAvg, twoDiffAvg, fourDiffAvg, targetLineLength);
    }
}

void vectorSpread(geoSolid& shape, int p1Index, int p2Index, double targetLineLength) 
{
    //get point coordinates of each index
    double* p1 = (shape.points + (3 * p1Index));
    double* p2 = (shape.points + (3 * p2Index));
    //average to find center
    double center[3] = {
        (p1[0] + p2[0]) / 2,
        (p1[1] + p2[1]) / 2,
        (p1[2] + p2[2]) / 2
    };
    //get distance from center to point (same for both)
    double twoPointDistance = distance(p1, center);
    //calculate vectors for each point
    double vectorCoefficient1[3] = {
        ((p1[0] - center[0])/twoPointDistance)*(targetLineLength/2-twoPointDistance),
        ((p1[1] - center[1])/twoPointDistance)*(targetLineLength/2-twoPointDistance),
        ((p1[2] - center[2])/twoPointDistance)*(targetLineLength/2-twoPointDistance)
    };
    double vectorCoefficient2[3] = {
        ((p2[0] - center[0]) / twoPointDistance)* (targetLineLength / 2 - twoPointDistance),
        ((p2[1] - center[1]) / twoPointDistance)* (targetLineLength / 2 - twoPointDistance),
        ((p2[2] - center[2]) / twoPointDistance)* (targetLineLength / 2 - twoPointDistance)
    };
    //apply vectors to each point
    double q1[] = { p1[0] + vectorCoefficient1[0],p1[1] + vectorCoefficient1[1],p1[2] + vectorCoefficient1[2] };
    double q2[] = { p2[0] + vectorCoefficient2[0],p2[1] + vectorCoefficient2[1],p2[2] + vectorCoefficient2[2] };
    for (int i = 0; i < 3;i++) {
        //apply adjusted points to shape points
        *(shape.points + (3 * p1Index) + i) = q1[i];
        *(shape.points + (3 * p2Index) + i) = q2[i];

    }

}

void fixRadialDistance(geoSolid& shape, double targetDistance) {
    //spread points out from the center as they have a tendency to collapse twoards 0 due to floating point math
    double center[3] = { 0,0,0 };
    for (int i = 0; i < shape.numPoints; i++) {
        double* point = (shape.points + (i * 3));
        double radialDistance = distance(point, center);
        double thetas[3];
        for (int u = 0; u<3;u++)
        {
            //I struggled to craft an equation that preserved the original sign of the coordinate so I just record the sign first, take the absolute value, and apply it back at the end
            int sign = 1;
            if (point[u] < 0) {
                sign *= -1;
            }
            //i'm not quite sure why this works, because its using the 3d distance and applying it to the 2d coordinates of the point, but it does the job
            *(shape.points + (i*3) + u) = std::abs(targetDistance * std::sin(std::asin(std::abs(point[u]) / radialDistance))) * sign;
        }
    }
}

void scale_points(geoSolid& shape, double scaler)
{ /* Scales up each value in the points of the 3d shape on startup*/
    for (int i = 0; i < shape.numPoints; i++) {
        for (int u = 0; u < 3; u++) {
            *(shape.points + (i * 3) + u) *= scaler;
        }
    }
}


/*######################## Helper Functions ########################*/
void arr_fill(double* target, double* src, int rows, int cols) {
    //fills in one array with the values of another
    for (int i = 0; i < rows; i++) {
        for (int u = 0; u < cols; u++) {
            *(target + (i * cols) + u) = *(src + (i * cols) + u);
        }
    }
}

void arr_fill(int* target, int* src, int rows, int cols) {
    //fills in one array with the values of another
    for (int i = 0; i < rows; i++) {
        for (int u = 0; u < cols; u++) {
            *(target + (i * cols) + u) = *(src + (i * cols) + u);
        }
    }
}

double distance(double* p1, double* p2) {
    //3d distance calculation
    return std::sqrt(
        std::pow(p1[0] - p2[0], 2) +
        std::pow(p1[1] - p2[1], 2) +
        std::pow(p1[2] - p2[2], 2));
}

void addFaceToPrint(geoSolid& shape, int face) {
    //add a face to shape.facesToPrint, checking to avoid duplicates
    bool faceAlreadyPresent = false;
    for (int i = 0; i < shape.numFacesToPrint; i++) {
        if (*(shape.facesToPrint + i) == face) faceAlreadyPresent = true;
    }
    if (!faceAlreadyPresent) {
        shape.facesToPrint[shape.numFacesToPrint] = face;
        shape.numFacesToPrint++;
    }
}

bool inBounds(double num, double bound1, double bound2) {
    //check if num is in the bounds bound1 & bound2
    if (num > bound1 && num < bound2) {
        return true;
    }
    if (num < bound1 && num > bound2) {
        return true;
    }
    return false;
}

int indexOfCharArr(char elem, char* arr, int len) {
    //get the index of char elem in char* arr
    for (int i = 0; i < len; i++) {
        if (*(arr + i) == elem) return i;
    }
    return -1;
}

int* indexOfIntArr(int* elem, int* arr, int rowlen, int collen, int elemlen) {
    //get the indecies of 2d int* arr where the ints from int*elem are all present in the element of arr, doesn't do duplicate checks but is unneccessary in this context 
    int* indecies = new int[2];
    int resultCount = 0;
    for (int i = 0; i < rowlen; i++) {
        bool matchedALl = true;
        for (int u = 0; u < elemlen; u++) {
            int lookingFor = *(elem + u);
            bool matchedThis = false;
            for (int w = 0; w < collen; w++) {
                if (*(arr + (i * collen) + w) == lookingFor) {
                    matchedThis = true;
                }
            }
            if (!matchedThis) {
                matchedALl = false;
            }
        }
        if (matchedALl) {
            indecies[resultCount] = i;
            resultCount++;
        }
    }
    return indecies;
}

bool checkConnectionGoodToPrint(geoSolid& shape, int connectionIndex) {
    //check if the connection is part of the faces in shape.facesToPrint
    int* faces = indexOfIntArr(shape.connectedPoints + (connectionIndex * 2), shape.faces, 6, 4, 2);
    for (int i = 0; i < shape.numFacesToPrint; i++) {
        if (faces[0] == *(shape.facesToPrint + i)) return true;
        if (faces[1] == *(shape.facesToPrint + i)) return true;
    }
    return false;
}

int* get_screen_dimensions(int maxWidth, int maxHeight) {
    //get the dimensions of the users screen
    int Yrange = maxHeight;
    int Xrange = maxWidth;
    int dimY = -1;
    int dimX = -1;
    std::cout << "\n\n\n\n\nPlease Maximize your window then press enter" << '\n';
    std::cin.ignore();
    //debug information forces default dimensions for ease of use while testing
    if (!print_debug_information) {
        std::string ymessage = "What is the highest number you can see on the screen? (0 for default values)";
        bool hasTriedY = false;
        bool hasTriedX = false;
        do {
            if (hasTriedY) { //if its not the first attempted input
                ymessage = "What is the highest number you can see on the screen? (must be an integer > 0)";
            }
            Yrange = maxHeight; //print a descending list of integers
            while (Yrange > 0) {
                std::cout << Yrange << '\n';
                Yrange--;
            }
            std::cout << ymessage << '\n';
            std::cin >> dimY; //get user input
            hasTriedY = true;
        } while (dimY < 0); //input validations
    }
    if (dimY == 0 || print_debug_information) { /*=======Default Values=======*/
        dimY = 61;
        dimX = 200;
    }
    else { /*=======Custom Values=======*/
        do { /* Input Validation Loop */
            //Display horizontal number line
            int Xcounter = 0;
            std::string xMeasureLines[3];
            do { //fills string array with integers written vertically
                std::string intString = std::to_string(Xcounter);
                int intStrLen = intString.length();
                for (int i = 0; i < 3; i++) {
                    xMeasureLines[i] += (i < intStrLen) ? intString[i] : ' ';
                }
                Xcounter++;
            } while (Xcounter <= Xrange);
            for (std::string line : xMeasureLines) { //prints strings
                std::cout << line << '\n';
            }

            //Get Input
            std::cout << "What is the right-most number you can see on the screen?" << '\n';
            std::cin >> dimX;
        } while (dimX < 0);
    }

    //return data
    int* return_pointer = new int[2];
    return_pointer[0] = dimX;
    return_pointer[1] = dimY + 1;
    return return_pointer;
}


/*######################## Display Functions ########################*/
void print_startup_banner() {
    std::string banner_lines[] = {
    "  _____      _",
    " /  __/_   _| |_   __       _ ",
    "|  |__| |_| |   \\/  __\\|\\ / /",
    " \\ ___\\_____|____/\\__/  \\  /",
    "                        |_/   ",
    "c++ final project by Justin Brown"

    };
    for (std::string line : banner_lines) {
        std::cout << line << '\n';
    }
}

void print(std::string* canvas, displayConfig& dispConf, int lowestY) {
    //create a single long string by concatinating the canvas strings, print them all out at once
    static int fixed_altitude = 0;
    std::string linestr;
    for (int i = 0; i < lowestY-dispConf.pushUp; i++) {
        linestr += '\n';
    }
    for (int i = dispConf.gridHeight - 1; i >= lowestY; i--) { //print canvas from top to bottom, i.e. canvas[0] should be last
        linestr += "\n" + canvas[i];
    }
    for (int i = 0; i < dispConf.pushUp; i++) {
        linestr += '\n';
    }
    std::cout << linestr;
}

void fill_newlines(int count) {
    //needed because before the whole vertical space of the terminal is filled there is a annoying jittering effect when things are printed quickly
    for (int i = 0; i < count;i++) {
        std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
    }
}

void print_points(double* pts) {
    //unused, just prints the point array, depricated and only here in case I want to use it again
    for (int c = 0; c < 8; c++) {
        std::cout << '\n';
        for (int i = 0; i < 2; i++) {
            if (i == 0) {
                std::cout << "{";
            }
            std::cout << *(pts + c * 2 + i);
            if (i == 2) {
                std::cout << "}";
            }
            else {
                std::cout << ", ";
            }
        }

    }
}


/*######################## Struct Initializers ########################*/
void init_shape(int shapeNum, geoSolid& shape) {
    //switch statement in case I want to add more shapes, unlikely
    switch (shapeNum) {
    case (0):
        shape.numPoints = 8;
        shape.numConnectedPoints = 12;
        shape.points = new double[8 * 3];
        double pts[][3] = {
        {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
        {-1, -1, 1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 1}
        };
        arr_fill(shape.points, *pts, 8, 3);
        int conpts[][2] = {
        {0, 1}, {0, 2}, {1, 3}, {2, 3},
        {0, 4}, {1, 5}, {2, 6}, {3, 7},
        {4, 5}, {4, 6}, {5, 7}, {6, 7}
        };
        shape.connectedPoints = new int[12 * 2];
        arr_fill(shape.connectedPoints, *conpts, 12, 2);
        int faces[][4] = { {0,1,2,3}, {0,1,4,5}, {0,2,4,6}, {1,5,3,7}, {4,5,6,7}, {2,3,6,7} };
        shape.numFaces = 6;
        shape.faces = new int[6 * 4];
        arr_fill(shape.faces, *faces, 6, 4);
        shape.projection = new int[3*8*3];
        shape.facesToPrint = new int[6];
        shape.lines = new double[12 * 4];
        
        break;
    }
}

void init_perspectives(perspectivePoint* perspectives) {
    //values for the sun & eye perspectivePoint instances, just a lot of text I wanted to take out of the main function
    perspectives[0].coord1 = 0;
    perspectives[0].coord2 = 2;
    perspectives[0].planeVal = -70;
    perspectives[0].planeIndex = 1;
    perspectives[0].points = new double[3];
    perspectives[0].points[0] = 25;
    perspectives[0].points[1] = -100;
    perspectives[0].points[2] = 0;

    perspectives[1].coord1 = 0;
    perspectives[1].coord2 = 1;
    perspectives[1].planeVal = 50;
    perspectives[1].planeIndex = 2;
    perspectives[1].points = new double[3];
    perspectives[1].points[0] = 0;
    perspectives[1].points[1] = 0;
    perspectives[1].points[2] = 200;
}
