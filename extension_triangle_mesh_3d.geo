SetFactory("OpenCASCADE");

Sphere(1) = {0, 0, 0, 0.5, 0, Pi/2, Pi/2};
Sphere(2) = {0, 0, 0, 1, 0, Pi/2, Pi/2};

BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };

Sphere(4) = {0, 0, 0, 0.5, 0, Pi/2, Pi/2};

Physical Volume(5) = {3,4};

Characteristic Length{ PointsOf{ Volume{3,4}; } } = 0.2;

