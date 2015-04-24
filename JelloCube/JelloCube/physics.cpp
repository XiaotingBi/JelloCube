/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <vector>
#include <iostream>

using namespace std;

double original_spring_length = 1.0 / 7.0;
// initially, the cube is on which side of the inclined plane.
int side;

// compute dot production of 2 vectors
double Vector_Dot_Product(struct point a, struct point b)
{
	double res = a.x * b.x + a.y * b.y + a.z * b.z;

	return res;
}

struct point externalForce(struct world * jello, int p_x, int p_y, int p_z)
{
	struct point ext_force;
	ext_force.x = 0;
	ext_force.y = 0;
	ext_force.z = 0;

	point this_point = jello->p[p_x][p_y][p_z];

	// the point is outside the bounding box
	if (this_point.x>=2 || this_point.x <= -2 || this_point.y >= 2 || this_point.y <= -2 || this_point.z >= 2 || this_point.z <= -2)
		return ext_force;

	double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0; // find out the force cube that contains this point
	// index in the array forceField
	int i, j, k;

	// x value
	for (i = 0; i < jello->resolution-1; i ++)
	{
		double xtemp0 = -2 + 4*(1.0 * i / (jello->resolution-1));
		double xtemp1 = -2 + 4*(1.0 * (i+1) / (jello->resolution-1));

		if (this_point.x >= xtemp0 && this_point.x <= xtemp1)
		{
			xmin = xtemp0;
			xmax = xtemp1;
			break;
		}
	}
	// z value
	for (k = 0; k < jello->resolution-1; k ++)
	{
		double ztemp0 = -2 + 4*(1.0 * k / (jello->resolution-1));
		double ztemp1 = -2 + 4*(1.0 * (k+1) / (jello->resolution-1));

		if (this_point.z >= ztemp0 && this_point.z <= ztemp1)
		{
			zmin = ztemp0;
			zmax = ztemp1;
			break;
		}
	}
	// y value
	for (j = 0; j < jello->resolution-1; j ++)
	{
		double ytemp0 = -2 + 4*(1.0 * j / (jello->resolution-1));
		double ytemp1 = -2 + 4*(1.0 * (j+1) / (jello->resolution-1));

		if (this_point.y >= ytemp0 && this_point.y <= ytemp1)
		{
			ymin = ytemp0;
			ymax = ytemp1;
			break;
		}
	}	

	// assign (x,y,z) values to 8 points of the cube
	struct point a, b, c, d, e, f, g, h;
	pMAKE(xmin, ymin, zmin, a);
	pMAKE(xmin, ymin, zmax, b);
	pMAKE(xmin, ymax, zmin, c);
	pMAKE(xmin, ymax, zmax, d);
	pMAKE(xmax, ymin, zmin, e);
	pMAKE(xmax, ymin, zmax, f);
	pMAKE(xmax, ymax, zmin, g);
	pMAKE(xmax, ymax, zmax, h);
	// interpolation
	int a_index = i*jello->resolution*jello->resolution+j*jello->resolution+k;
	int b_index = i*jello->resolution*jello->resolution+j*jello->resolution+k+1;
	int c_index = i*jello->resolution*jello->resolution+(j+1)*jello->resolution+k;
	int d_index = i*jello->resolution*jello->resolution+(j+1)*jello->resolution+k+1;
	int e_index = (i+1)*jello->resolution*jello->resolution+j*jello->resolution+k;
	int f_index = (i+1)*jello->resolution*jello->resolution+j*jello->resolution+k+1;
	int g_index = (i+1)*jello->resolution*jello->resolution+(j+1)*jello->resolution+k;
	int h_index = (i+1)*jello->resolution*jello->resolution+(j+1)*jello->resolution+k+1;

	// force field value
	point a_force = jello->forceField[a_index];
	point b_force = jello->forceField[b_index];
	point c_force = jello->forceField[c_index];
	point d_force = jello->forceField[d_index];
	point e_force = jello->forceField[e_index];
	point f_force = jello->forceField[f_index];
	point g_force = jello->forceField[g_index];
	point h_force = jello->forceField[h_index];
	
	double x, y, z;
	x = abs((this_point.x - xmin) / (xmax - xmin));
	y = abs((this_point.y - ymin) / (ymax - ymin));
	z = abs((this_point.z - zmin) / (zmax - zmin));

	ext_force.x = (1-x)*(1-y)*(1-z)*a_force.x + (1-x)*(1-y)*z*b_force.x 
				+ (1-x)*y*(1-z)*c_force.x + (1-x)*y*z*d_force.x
				+ x*(1-y)*(1-z)*e_force.x + x*(1-y)*z*f_force.x
				+ x*y*(1-z)*g_force.x + x*y*z*h_force.x;
	ext_force.y = (1-x)*(1-y)*(1-z)*a_force.y + (1-x)*(1-y)*z*b_force.y
				+ (1-x)*y*(1-z)*c_force.y + (1-x)*y*z*d_force.y
				+ x*(1-y)*(1-z)*e_force.y + x*(1-y)*z*f_force.y
				+ x*y*(1-z)*g_force.y + x*y*z*h_force.y;
	ext_force.z = (1-x)*(1-y)*(1-z)*a_force.z + (1-x)*(1-y)*z*b_force.z
				+ (1-x)*y*(1-z)*c_force.z + (1-x)*y*z*d_force.z
				+ x*(1-y)*(1-z)*e_force.z + x*(1-y)*z*f_force.z
				+ x*y*(1-z)*g_force.z + x*y*z*h_force.z;

	return ext_force;
}

struct point externalForce2(struct world * jello, int p_x, int p_y, int p_z)
{
	struct point ext_force;
	ext_force.x = 0;
	ext_force.y = 0;
	ext_force.z = 0;

	point this_point = jello->p[p_x][p_y][p_z];

	// the point is outside the bounding box
	if (this_point.x>=2 || this_point.x <= -2 || this_point.y >= 2 || this_point.y <= -2 || this_point.z >= 2 || this_point.z <= -2)
		return ext_force;

	double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0; // find out the force cube that contains this point
	// index in the array forceField
	int i, j, k;

	// x value
	for (i = 0; i < jello->resolution-1; i ++)
	{
		double xtemp0 = -2 + 4*(1.0 * i / (jello->resolution-1));
		double xtemp1 = -2 + 4*(1.0 * (i+1) / (jello->resolution-1));

		if (this_point.x >= xtemp0 && this_point.x <= xtemp1)
		{
			xmin = xtemp0;
			xmax = xtemp1;
			break;
		}
	}
	// z value
	for (k = 0; k < jello->resolution-1; k ++)
	{
		double ztemp0 = -2 + 4*(1.0 * k / (jello->resolution-1));
		double ztemp1 = -2 + 4*(1.0 * (k+1) / (jello->resolution-1));

		if (this_point.z >= ztemp0 && this_point.z <= ztemp1)
		{
			zmin = ztemp0;
			zmax = ztemp1;
			break;
		}
	}
	// y value
	for (j = 0; j < jello->resolution-1; j ++)
	{
		double ytemp0 = -2 + 4*(1.0 * j / (jello->resolution-1));
		double ytemp1 = -2 + 4*(1.0 * (j+1) / (jello->resolution-1));

		if (this_point.y >= ytemp0 && this_point.y <= ytemp1)
		{
			ymin = ytemp0;
			ymax = ytemp1;
			break;
		}
	}	

	// assign (x,y,z) values to 8 points of the cube
	struct point a, b, c, d, e, f, g, h;
	pMAKE(xmin, ymin, zmin, a);
	pMAKE(xmin, ymin, zmax, b);
	pMAKE(xmin, ymax, zmin, c);
	pMAKE(xmin, ymax, zmax, d);
	pMAKE(xmax, ymin, zmin, e);
	pMAKE(xmax, ymin, zmax, f);
	pMAKE(xmax, ymax, zmin, g);
	pMAKE(xmax, ymax, zmax, h);
	// interpolation
	int a_index = i*jello->resolution*jello->resolution+j*jello->resolution+k;
	int b_index = i*jello->resolution*jello->resolution+j*jello->resolution+k+1;
	int d_index = i*jello->resolution*jello->resolution+(j+1)*jello->resolution+k;
	int c_index = i*jello->resolution*jello->resolution+(j+1)*jello->resolution+k+1;
	int e_index = (i+1)*jello->resolution*jello->resolution+j*jello->resolution+k;
	int f_index = (i+1)*jello->resolution*jello->resolution+j*jello->resolution+k+1;
	int h_index = (i+1)*jello->resolution*jello->resolution+(j+1)*jello->resolution+k;
	int g_index = (i+1)*jello->resolution*jello->resolution+(j+1)*jello->resolution+k+1;

	// force field value
	point a_force = jello->forceField[a_index];
	point b_force = jello->forceField[b_index];
	point c_force = jello->forceField[c_index];
	point d_force = jello->forceField[d_index];
	point e_force = jello->forceField[e_index];
	point f_force = jello->forceField[f_index];
	point g_force = jello->forceField[g_index];
	point h_force = jello->forceField[h_index];
	// point b, c
	double b_effect = 1 - abs((this_point.y - b.y) / (c.y - b.y));
	double c_effect = 1 - b_effect;
	point bc;
	pMAKE(b_effect*b_force.x+c_effect*c_force.x, b_effect*b_force.y+c_effect*c_force.y, b_effect*b_force.z+c_effect*c_force.z, bc);
	// point a, d
	double a_effect = 1 - abs((this_point.y - a.y) / (d.y - a.y));
	double d_effect = 1 - a_effect;
	point ad;
	pMAKE(a_effect*a_force.x+d_effect*d_force.x, a_effect*a_force.y+d_effect*d_force.y, a_effect*a_force.z+d_effect*d_force.z, ad);
	// point f, g
	double f_effect = 1 - abs((this_point.y - f.y) / (g.y - f.y));
	double g_effect = 1 - f_effect;
	point fg;
	pMAKE(f_effect*f_force.x+g_effect*g_force.x, f_effect*f_force.y+g_effect*g_force.y, f_effect*f_force.z+g_effect*g_force.z, fg);
	// point e, h
	double e_effect = 1 - abs((this_point.y - e.y) / (h.y - e.y));
	double h_effect = 1 - e_effect;
	point eh;
	pMAKE(e_effect*e_force.x+h_effect*h_force.x, e_effect*e_force.y+h_effect*h_force.y, e_effect*e_force.z+h_effect*h_force.z, eh);
	// face abcd
	double abcd_bc = 1 - abs((b.z - this_point.z) / (a.z - b.z));
	double abcd_ad = 1 - abcd_bc;
	point abcd_face;
	pMAKE(bc.x*abcd_bc+ad.x*abcd_ad, bc.y*abcd_bc+ad.y*abcd_ad, bc.z*abcd_bc+ad.z*abcd_ad, abcd_face);
	// face efgh
	double efgh_fg = 1 - abs((f.z - this_point.z) / (e.z - f.z));
	double efgh_eh = 1 - efgh_fg;
	point efgh_face;
	pMAKE(fg.x*efgh_fg+efgh_eh*eh.x, fg.y*efgh_fg+eh.y*efgh_eh, fg.z*efgh_fg+eh.z*efgh_eh, efgh_face);

	double final_abcd = 1 - abs((this_point.x - a.x) / (a.x - e.x));
	double final_efgh = 1 - final_abcd;
	ext_force.x = final_abcd*abcd_face.x + final_efgh*efgh_face.x;
	ext_force.y = final_abcd*abcd_face.y + final_efgh*efgh_face.y;
	ext_force.z = final_abcd*abcd_face.z + final_efgh*efgh_face.z;

	/*cout << ext_force.x << ", " << ext_force.y << ", " << ext_force.z << endl;*/

	return ext_force;
}


struct point bendForce(struct world * jello, int p_x, int p_y, int p_z)
{
	struct point bend_force;
	bend_force.x = 0;
	bend_force.y = 0;
	bend_force.z = 0;

	vector<struct point> points_index; // use struct point to save the index (i, j, k) of a point

	struct point temp;
	temp.x = 0;
	temp.y = 0;
	temp.z = 0;
	if (p_x > 1)
	{	
		temp.x = p_x - 2;
		temp.y = p_y;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_x < 6)
	{
		temp.x = p_x + 2;
		temp.y = p_y;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_y > 1)
	{
		temp.x = p_x;
		temp.y = p_y - 2;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_y < 6)
	{
		temp.x = p_x;
		temp.y = p_y + 2;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_z > 1)
	{
		temp.x = p_x;
		temp.y = p_y;
		temp.z = p_z - 2;
		points_index.push_back(temp);
	}
	if (p_z < 6)
	{
		temp.x = p_x;
		temp.y = p_y;
		temp.z = p_z + 2;
		points_index.push_back(temp);
	} 

	struct point sum_force; // sum of all kinds of forces
	sum_force.x = 0;
	sum_force.y = 0;
	sum_force.z = 0;
	// calculate forces imposed on (p_i, p_j, p_k) by its neighbors
	for (int i = 0; i < points_index.size(); i ++)
	{
		// (x, y, z) of a point
		int i_x = points_index[i].x;
		int i_y = points_index[i].y;
		int i_z = points_index[i].z;

		struct point L;
		pDIFFERENCE(jello->p[p_x][p_y][p_z], jello->p[i_x][i_y][i_z], L); // position vector: p(i,j,k)-p(ii,jj,kk)
		struct point diff_of_vel;
		pDIFFERENCE(jello->v[p_x][p_y][p_z], jello->v[i_x][i_y][i_z], diff_of_vel); // velocity vector: v(i,j,k)-v(ii,jj,kk)

		double length_of_L = sqrt(L.x * L.x + L.y * L.y + L.z * L.z); // length of L
		struct point hook_force; // Hook Force
		pMULTIPLY(L, -1*jello->kElastic*(length_of_L-original_spring_length*2)/length_of_L, hook_force);
		pSUM(sum_force, hook_force, sum_force);

		struct point damping_force; // Damping Force
		double temp_damp = Vector_Dot_Product(L, diff_of_vel); // L * (v(i,j,k)-v(ii,jj,kk))
		pMULTIPLY(L, -1*jello->dElastic*temp_damp/length_of_L/length_of_L, damping_force);
		pSUM(sum_force, damping_force, sum_force);
	}
	bend_force.x = sum_force.x;
	bend_force.y = sum_force.y;
	bend_force.z = sum_force.z;

	return bend_force;
}

struct point shearForce(struct world * jello, int p_x, int p_y, int p_z)
{
	struct point shear_force;
	shear_force.x = 0;
	shear_force.y = 0;
	shear_force.z = 0;

	vector<struct point> points_index; // use struct point to save the index (i, j, k) of a point
	
	struct point temp;
	temp.x = 0;
	temp.y = 0;
	temp.z = 0;

	// same x value

	// x, y, z-1: this situation does not exist
	// x, y, z+1: this situation does not exist
	if (p_y > 0)
	{
		// same z value
		// x, y-1, z: this situation does not exist
		// x, y-1, z-1
		if (p_z > 0)
		{
			temp.x = p_x;
			temp.y = p_y - 1;
			temp.z = p_z - 1;
			points_index.push_back(temp);
		}
		// x, y-1, z+1
		if (p_z < 7)
		{
			temp.x = p_x;
			temp.y = p_y - 1;
			temp.z = p_z + 1;
			points_index.push_back(temp);
		}		
	}
	if (p_y < 7)
	{
		// same z value
		// x, y+1, z
		// x, y+1, z-1
		if (p_z > 0)
		{
			temp.x = p_x;
			temp.y = p_y + 1;
			temp.z = p_z - 1;
			points_index.push_back(temp);
		}
		// x, y+1, z+1
		if (p_z < 7)
		{
			temp.x = p_x;
			temp.y = p_y + 1;
			temp.z = p_z + 1;
			points_index.push_back(temp);
		}	
	}
	
	if (p_x > 0)
	{
		// same y value: y = 0 || y = 7
		// x-1, y, z
		// x-1, y, z-1
		if (p_z > 0)
		{	
			temp.x = p_x - 1;
			temp.y = p_y;
			temp.z = p_z - 1;
			points_index.push_back(temp);
		}
		// x-1, y, z+1
		if (p_z < 7)
		{	
			temp.x = p_x - 1;
			temp.y = p_y;
			temp.z = p_z + 1;
			points_index.push_back(temp);
		}

		if (p_y > 0)
		{
			// same z value
			// x-1, y-1, z
			temp.x = p_x - 1;
			temp.y = p_y - 1;
			temp.z = p_z;
			points_index.push_back(temp);
			// x-1, y-1, z-1
			if (p_z > 0)
			{
				temp.x = p_x - 1;
				temp.y = p_y - 1;
				temp.z = p_z - 1;
				points_index.push_back(temp);
			}
			// x-1, y-1, z+1
			if (p_z < 7)
			{
				temp.x = p_x - 1;
				temp.y = p_y - 1;
				temp.z = p_z + 1;
				points_index.push_back(temp);
			}		
		}
		if (p_y < 7)
		{
			// same z value
			// x-1, y+1, z
			temp.x = p_x - 1;
			temp.y = p_y + 1;
			temp.z = p_z;
			points_index.push_back(temp);
			// x-1, y+1, z-1
			if (p_z > 0)
			{
				temp.x = p_x - 1;
				temp.y = p_y + 1;
				temp.z = p_z - 1;
				points_index.push_back(temp);
			}
			// x-1, y+1, z+1
			if (p_z < 7)
			{
				temp.x = p_x - 1;
				temp.y = p_y + 1;
				temp.z = p_z + 1;
				points_index.push_back(temp);
			}	
		}
	} 

	if (p_x < 7)
	{
		// same y value
		// x+1, y, z
		// x+1, y, z-1
		if (p_z > 0)
		{	
			temp.x = p_x + 1;
			temp.y = p_y;
			temp.z = p_z - 1;
			points_index.push_back(temp);
		}
		// x+1, y, z+1
		if (p_z < 7)
		{	
			temp.x = p_x + 1;
			temp.y = p_y;
			temp.z = p_z + 1;
			points_index.push_back(temp);
		}

		if (p_y > 0)
		{
			// same z value
			// x+1, y-1, z
			temp.x = p_x + 1;
			temp.y = p_y - 1;
			temp.z = p_z;
			points_index.push_back(temp);
			// x+1, y-1, z-1
			if (p_z > 0)
			{
				temp.x = p_x + 1;
				temp.y = p_y - 1;
				temp.z = p_z - 1;
				points_index.push_back(temp);
			}
			// x+1, y-1, z+1
			if (p_z < 7)
			{
				temp.x = p_x + 1;
				temp.y = p_y - 1;
				temp.z = p_z + 1;
				points_index.push_back(temp);
			}		
		}
		if (p_y < 7)
		{
			// same z value
			// x+1, y+1, z
			temp.x = p_x + 1;
			temp.y = p_y + 1;
			temp.z = p_z;
			points_index.push_back(temp);
			// x+1, y+1, z-1
			if (p_z > 0)
			{
				temp.x = p_x + 1;
				temp.y = p_y + 1;
				temp.z = p_z - 1;
				points_index.push_back(temp);
			}
			// x+1, y+1, z+1
			if (p_z < 7)
			{
				temp.x = p_x + 1;
				temp.y = p_y + 1;
				temp.z = p_z + 1;
				points_index.push_back(temp);
			}	
		}
	} 
	
	struct point sum_force; // sum of all kinds of forces
	sum_force.x = 0;
	sum_force.y = 0;
	sum_force.z = 0;

	// plane diagonal spring
	double scalar = sqrt(2.0);
	// calculate forces imposed on (p_i, p_j, p_k) by its neighbors
	for (int i = 0; i < points_index.size(); i ++)
	{		
		// (x, y, z) of a point
		int i_x = points_index[i].x;
		int i_y = points_index[i].y;
		int i_z = points_index[i].z;

		// cube diagonal spring
		if (abs(p_y - i_y) == 1 && abs(p_z - i_z) == 1 && abs(p_x - i_x) == 1)
			scalar = sqrt(3.0);
		else
			scalar = sqrt(2.0);

		struct point L;
		pDIFFERENCE(jello->p[p_x][p_y][p_z], jello->p[i_x][i_y][i_z], L); // position vector: p(i,j,k)-p(ii,jj,kk)
		struct point diff_of_vel;
		pDIFFERENCE(jello->v[p_x][p_y][p_z], jello->v[i_x][i_y][i_z], diff_of_vel); // velocity vector: v(i,j,k)-v(ii,jj,kk)

		double length_of_L = sqrt(L.x * L.x + L.y * L.y + L.z * L.z); // length of L
		struct point hook_force; // Hook Force
		pMULTIPLY(L, -1*jello->kElastic*(length_of_L-original_spring_length*scalar)/length_of_L, hook_force);
		pSUM(sum_force, hook_force, sum_force);

		struct point damping_force; // Damping Force
		double temp_damp = Vector_Dot_Product(L, diff_of_vel); // L * (v(i,j,k)-v(ii,jj,kk))
		pMULTIPLY(L, -1*jello->dElastic*temp_damp/length_of_L/length_of_L, damping_force);
		pSUM(sum_force, damping_force, sum_force);
	}
	shear_force.x = sum_force.x;
	shear_force.y = sum_force.y;
	shear_force.z = sum_force.z;

	return shear_force;
}

struct point structForce(struct world * jello, int p_x, int p_y, int p_z)
{	
	struct point struct_force;
	struct_force.x = 0;
	struct_force.y = 0;
	struct_force.z = 0;

	vector<struct point> points_index; // use struct point to save the index (i, j, k) of a point

	struct point temp;
	temp.x = 0;
	temp.y = 0;
	temp.z = 0;
	if (p_x > 0)
	{	
		temp.x = p_x - 1;
		temp.y = p_y;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_x < 7)
	{
		temp.x = p_x + 1;
		temp.y = p_y;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_y > 0)
	{
		temp.x = p_x;
		temp.y = p_y - 1;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_y < 7)
	{
		temp.x = p_x;
		temp.y = p_y + 1;
		temp.z = p_z;
		points_index.push_back(temp);
	}
	if (p_z > 0)
	{
		temp.x = p_x;
		temp.y = p_y;
		temp.z = p_z - 1;
		points_index.push_back(temp);
	}
	if (p_z < 7)
	{
		temp.x = p_x;
		temp.y = p_y;
		temp.z = p_z + 1;
		points_index.push_back(temp);
	} 

	struct point sum_force; // sum of all kinds of forces
	sum_force.x = 0;
	sum_force.y = 0;
	sum_force.z = 0;
	// calculate forces imposed on (p_i, p_j, p_k) by its neighbors
	for (int i = 0; i < points_index.size(); i ++)
	{
		// (x, y, z) of a point
		int i_x = points_index[i].x;
		int i_y = points_index[i].y;
		int i_z = points_index[i].z;

		struct point L;
		pDIFFERENCE(jello->p[p_x][p_y][p_z], jello->p[i_x][i_y][i_z], L); // position vector: p(i,j,k)-p(ii,jj,kk)
		struct point diff_of_vel;
		pDIFFERENCE(jello->v[p_x][p_y][p_z], jello->v[i_x][i_y][i_z], diff_of_vel); // velocity vector: v(i,j,k)-v(ii,jj,kk)

		double length_of_L = sqrt(L.x * L.x + L.y * L.y + L.z * L.z); // length of L
		struct point hook_force; // Hook Force
		pMULTIPLY(L, -1*jello->kElastic*(length_of_L-original_spring_length)/length_of_L, hook_force);
		pSUM(sum_force, hook_force, sum_force);
		
		struct point damping_force; // Damping Force
		double temp_damp = Vector_Dot_Product(L, diff_of_vel); // L * (v(i,j,k)-v(ii,jj,kk))
		pMULTIPLY(L, -1*jello->dElastic*temp_damp/length_of_L/length_of_L, damping_force);
		pSUM(sum_force, damping_force, sum_force);
	}
	struct_force.x = sum_force.x;
	struct_force.y = sum_force.y;
	struct_force.z = sum_force.z;

	return struct_force;
}

struct point collideForce(struct world * jello, int p_x, int p_y, int p_z)
{
	point collision_force;
	collision_force.x = 0;
	collision_force.y = 0;
	collision_force.z = 0;
	// the real position which may be outside the bounding box; if not, no collision force
	point real_position = jello->p[p_x][p_y][p_z];

	point bounce_direction; // normal vector
	bounce_direction.x = 0;
	bounce_direction.y = 0;
	bounce_direction.z = 0;
	/* the amount of penetration to x,y,z axes */
	double extra_x = 0, extra_y = 0, extra_z = 0;
	/* velocity when colliding */
	point collide_velocity = jello->v[p_x][p_y][p_z];

	// collision detection, to detect the box is outside of which faces (may be at most 3 faces)
	if ((real_position.x >= 2) || (real_position.x < -2))
	{		
		if (real_position.x >= 2)
		{		
			// magnitude of the force
			extra_x = real_position.x - 2;
			// direction of the bouncing off force
			bounce_direction.x = -1;
			bounce_direction.y = 0;
			bounce_direction.z = 0;		
		} else 
		{
			extra_x = -2 - real_position.x;
			/*extra_x *= -1;*/
			// direction of the bouncing off force
			bounce_direction.x = 1;
			bounce_direction.y = 0;
			bounce_direction.z = 0;	
		}
		// Hook Law used in collision
		point hook_force_collision;
		double t_hook = jello->kCollision * extra_x;
		pMULTIPLY(bounce_direction, t_hook, hook_force_collision);
		pSUM(collision_force, hook_force_collision, collision_force);

		// Damping force in collision
		point damp_force_colllision;
		double t_damp = Vector_Dot_Product(bounce_direction, collide_velocity);
		t_damp = t_damp * jello->dCollision;
		pMULTIPLY(bounce_direction, t_damp, damp_force_colllision);
		pSUM(collision_force, damp_force_colllision, collision_force);
	}

	if ((real_position.y >= 2) || (real_position.y <= -2))
	{		
		if (real_position.y >= 2)
		{	// magnitude of the force
			extra_y = real_position.y - 2;
			// direction of the bouncing off force
			bounce_direction.x = 0;
			bounce_direction.y = -1;
			bounce_direction.z = 0;
		} else
		{
			extra_y = -2 - real_position.y;
			bounce_direction.x = 0;
			bounce_direction.y = 1;
			bounce_direction.z = 0;
		}
		// Hook Law used in collision
		point hook_force_collision;
		double t_hook = jello->kCollision * extra_y;
		pMULTIPLY(bounce_direction, t_hook, hook_force_collision);
		pSUM(collision_force, hook_force_collision, collision_force);

		// Damping force in collision
		point damp_force_colllision;
		double t_damp = Vector_Dot_Product(bounce_direction, collide_velocity);
		t_damp = t_damp * jello->dCollision;
		pMULTIPLY(bounce_direction, t_damp, damp_force_colllision);
		pSUM(collision_force, damp_force_colllision, collision_force);
	} 

	if ((real_position.z >= 2) || (real_position.z <= -2))
	{
		if (real_position.z >= 2)
		{
			// magnitude of the force
			extra_z = real_position.z - 2;
			bounce_direction.x = 0;
			bounce_direction.y = 0;
			bounce_direction.z = -1;
		} else
		{
			// magnitude of the force
			extra_z = -1 * real_position.z - 2;
			bounce_direction.x = 0;
			bounce_direction.y = 0;
			bounce_direction.z = 1;
		}
		// Hook Law used in collision
		point hook_force_collision;
		double t_hook = jello->kCollision * extra_z;
		pMULTIPLY(bounce_direction, t_hook, hook_force_collision);
		pSUM(collision_force, hook_force_collision, collision_force);

		// Damping force in collision
		point damp_force_colllision;
		double t_damp = Vector_Dot_Product(bounce_direction, collide_velocity);
		t_damp = t_damp * jello->dCollision;
		pMULTIPLY(bounce_direction, t_damp, damp_force_colllision);
		pSUM(collision_force, damp_force_colllision, collision_force);
	}

	return collision_force;
}
/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
	for (int i = 0; i < 8; i ++)
		for (int j = 0; j < 8; j ++)
			for (int k = 0; k < 8; k ++)
			{				
				point struct_force = structForce(jello, i, j, k);
				point shear_force = shearForce(jello, i, j, k);
				point bend_force = bendForce(jello, i, j, k);				
				point collide_force = collideForce(jello, i, j, k);				
				point external_force = externalForce(jello, i, j, k); 

				point total_force;
 				pSUM(struct_force, bend_force, total_force);
 				pSUM(total_force, shear_force, total_force);
				pSUM(total_force, collide_force, total_force);
				pSUM(total_force, external_force, total_force);			

				pMULTIPLY(total_force, 1.0 / jello->mass, a[i][j][k]);
			}
			
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}