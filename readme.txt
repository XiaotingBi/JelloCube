Features:

1. Animate the movement of a jello cube based on a realistic physical model.
	There are three kinds of springs: structural springs, shear springs and bend springs.
	Internal force of a spring: Hook's force and damping force.
2. Incorporate a bounding box, including collision detection and response.
	When a point is outside the bounding box, there will be a penalty force to push the point back into the bounding box.
	The collision force also contains hook's force and damping force.
3. Implement an arbitrary non-homogeneous time-independent force field, with appropriate force-field interpolation.
	The formula to interpolate force field:
	F(x,y,z) = (1-x)(1-y)(1-z)F(0,0,0) + (1-x)(1-y)zF(0,0,1)
			 + (1-x)y(1-z)F(0,1,0) + (1-x)yzF(0,1,1)
			 + x(1-y)(1-z)F(1,0,0) + x(1-y)zF(1,0,1)
			 + xy(1-z)F(1,1,0) + xyzF(1,1,1);
4. Use the provided code to interactively position the camera, use proper lighting, and render the cube in wireframe and triangle mode.
5. The program works with arbitrary world files.