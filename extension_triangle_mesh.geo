//Merge "simplex_body_fitting0.iges";
//Line Loop (1) = {1};
//Characteristic Length { 1 } = 0.05;
//Plane Surface(1) = {1};
Point (2) = {0, 0, 0 ,0.1 };
Point (3) = {2, 0, 0 ,0.1 };
Point (4) = {0, 2, 0 ,0.1 };
Point (5) = {2, 2, 0 ,0.1 };
Point (6) = {1, 0, 0 ,0.1 };
Point (7) = {0, 1, 0 ,0.1 };

Line (2) = {2, 7 };
Line (3) = {7, 4 };
Line (4) = {4, 5 };
Line (5) = {5, 3 };
Line (6) = {3, 6 };
Line (7) = {6, 2 };

//+
Circle(8) = {6, 2, 7};
//+
Curve Loop(1) = {7, 2, -8};
//+
Plane Surface(1) = {1};


Curve Loop(2) = {2, 3, 4, 5, 6, 7}; 
Physical Line(1) = {2, 3, 4, 5, 6, 7};
Plane Surface(2 ) = {1, 2 };
Plane Surface(3 ) = {1};
Physical Surface(4 ) = { 2,3 };

