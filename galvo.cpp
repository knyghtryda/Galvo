#include "galvo.h"
#include "arduino.h"
#include "mult16x8.h"
#include "mult16x16.h"

Galvo::Galvo() {
	//Setting some sane values
	setMMSize(150, 150, 150, 10);
	setScale(1, 1);
	Initialize();
}

Galvo::Galvo(float x_length, float y_length, float z_length, float e_length, float x_scale, float y_scale) {
	setMMSize(x_length, y_length, z_length, e_length);
	setScale(x_scale, y_scale);
	Initialize();
}

void Galvo::Initialize() {
	calcMinStepSize();
	calcMaxStepsPerUnit();
	calcMin();
	calcMax();
	calcT0Max();
	calcSize();
	calcCalStepSize();
	calcStepSize();
	calcTMax();
	calcStepsPerMM();
	calcZSize();
	e = mm_size[E] * steps_per_mm[X];
	position[X] = minVal[X];
	position[Y] = minVal[Y];
}

void Galvo::CalcCalibrationTable() {
	float x_pos = (float)minVal[X], y_pos = (float)minVal[Y];
	float x_tmp, y_tmp;
	for (int j = 0; j < points; j++) {
		for (int i = 0; i < points; i++) {
			x_tmp = x_pos - (float)center[X];
			y_tmp = y_pos - (float)center[Y];
			y_tmp = atan2(y_tmp, z_size[Y]) * t_max[Y]; // compute y first as its designated the independent (second galvo) axis
			x_tmp = atan2(x_tmp, sqrt((y_tmp * y_tmp) + ((float)z_size[X] * (float)z_size[X])) + (float)e) * t_max[X];  // compute x next as its dependent on the Y movement
			x_tmp += (float)center[X] - x_pos;
			y_tmp += (float)center[Y] - y_pos;
			offsets[i][j][X] = (int)x_tmp;
			offsets[i][j][Y] = (int)y_tmp;
			//printPair(offsets[i][j]);
			x_pos += cal_step_size[X];
			if (x_pos > maxVal[X]) x_pos = minVal[X];
		}
		y_pos += cal_step_size[Y];
		if (y_pos > maxVal[Y]) y_pos = minVal[Y];
		//Serial.println("");
	}
}

#ifdef TRIANGLE_MESH
/*
*	2x2 Mesh Grid Example
*	----------------------
*	All larger grids extrapolate directory from this example, and must have an
*	even number of grid spaces.
*
*	This is the orientation of the quadrants.  This may need to be standardized.
*		 |
*	  Q2 | Q1
*	-----------
*	  Q3 | Q4
*		 |
*
*	This is how the triangles are arranged in order to for the comparison math 
*	to be Straight forward
* (0,2)-----(1,2)-----(2,2)
*	| \	  S1  |	 S1	  /	|
*	|	 \	  |	   /	|
*	| S0	\ | /	 S0	|
* (0,1)-----(1,1)-----(2,1)
*	| S0	/ | \	 S0	|
*	|	 /	  |	   \	|
*	| /   S1  |	 S1	  \	|
* (0,0)-----(1,0)-----(2,0)
*
*	Q1 example 
*	01-------11
*	| S1	/ |
*	|	 /	  | <--- cal_step_size[Y]
*	| /   S0  |
*	00-------10
*		 ^
*		 |
*	cal_step_size[X]
*
*	S0 = abs(dx * css[Y]) > abs(dy * css[X]) 
*	S1 = abs(dx * css[Y]) <= abs(dy * css[X]) 
*	These multiplys are necessary in order to account for when cal_grid_size
*	is not equal between X and Y.
*	This applies for all quadrants, just rotated.
*/
void Galvo::CalcSlopeTable() {
	int s1[2][2], s2[2][2];
	int offset00[2], offset01[2], offset10[2], offset11[2];
	unsigned char quadrant;
	for (int j = 0; j < steps; j++) {
		for (int i = 0; i < steps; i++) {
			offset00[X] = offsets[i][j][X];
			offset01[X] = offsets[i][j + 1][X];
			offset10[X] = offsets[i + 1][j][X];
			offset11[X] = offsets[i + 1][j + 1][X];
			offset00[Y] = offsets[i][j][Y];
			offset01[Y] = offsets[i][j + 1][Y];
			offset10[Y] = offsets[i + 1][j][Y];
			offset11[Y] = offsets[i + 1][j + 1][Y];
			if (i < steps >> 1) {
				if (j < steps >> 1) {
					slopes[i][j][S0][X][X] = (int)(((float)(offset11[X] - offset01[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][X][Y] = (int)(((float)(offset01[X] - offset00[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][X] = (int)(((float)(offset11[Y] - offset01[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][Y] = (int)(((float)(offset01[Y] - offset00[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][X] = (int)(((float)(offset10[X] - offset00[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][Y] = (int)(((float)(offset11[X] - offset10[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][X] = (int)(((float)(offset10[Y] - offset00[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][Y] = (int)(((float)(offset11[Y] - offset10[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
				}
				else {
					slopes[i][j][S0][X][X] = (int)(((float)(offset10[X] - offset00[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][X][Y] = (int)(((float)(offset01[X] - offset00[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][X] = (int)(((float)(offset10[Y] - offset00[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][Y] = (int)(((float)(offset01[Y] - offset00[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][X] = (int)(((float)(offset11[X] - offset01[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][Y] = (int)(((float)(offset11[X] - offset10[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][X] = (int)(((float)(offset11[Y] - offset01[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][Y] = (int)(((float)(offset11[Y] - offset10[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
				}
			}
			else {
				if (j < steps >> 1) {
					slopes[i][j][S0][X][X] = (int)(((float)(offset11[X] - offset01[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][X][Y] = (int)(((float)(offset11[X] - offset10[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][X] = (int)(((float)(offset11[Y] - offset01[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][Y] = (int)(((float)(offset11[Y] - offset10[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][X] = (int)(((float)(offset10[X] - offset00[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][Y] = (int)(((float)(offset01[X] - offset00[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][X] = (int)(((float)(offset10[Y] - offset00[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][Y] = (int)(((float)(offset01[Y] - offset00[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
				}
				else {
					slopes[i][j][S0][X][X] = (int)(((float)(offset10[X] - offset00[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][X][Y] = (int)(((float)(offset11[X] - offset01[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][X] = (int)(((float)(offset10[Y] - offset00[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S0][Y][Y] = (int)(((float)(offset11[Y] - offset01[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][X] = (int)(((float)(offset11[X] - offset01[X]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][X][Y] = (int)(((float)(offset01[X] - offset00[X]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][X] = (int)(((float)(offset11[Y] - offset01[Y]) / (float)cal_step_size[X]) * (float)FIXED_DIV);
					slopes[i][j][S1][Y][Y] = (int)(((float)(offset01[Y] - offset00[Y]) / (float)cal_step_size[Y]) * (float)FIXED_DIV);
				}
			}
			/*
			printPair(offset01);
			printPair(offset11);
			Serial.println("");
			printPair(offset00);
			printPair(offset10);
			Serial.println("");
			Serial.print("Fixed point S0 vectors = X: ");
			printPair(slopes[i][j][S0][X]);
			Serial.print(", Y: ");
			printPair(slopes[i][j][S0][Y]);
			Serial.println("");
			Serial.print("Fixed point S1 vectors = X: ");
			printPair(slopes[i][j][S1][X]);
			Serial.print(", Y: ");
			printPair(slopes[i][j][S1][Y]);
			Serial.println("");
			*/
		}
	}
}

void Galvo::ResetSlopeTable() {
	memset(slopes, 0, sizeof(slopes));
}
#endif

void Galvo::ResetCalibrationTable() {
	memset(offsets, 0, sizeof(offsets));
}

//Applies the offset table to a set of coordinates
void Galvo::ApplyOffsets(unsigned int * val) {
	// shifting all coordinate values by 4 (divide by 16) in order to stop overflowing 32 bit longs
	// This reduces bilinear interpolation resolution, but is necessary in order to not use floats.  
	// This is not perfect and will still overflow if low values are chosen for number of grid points.

	unsigned int sx = cal_step_size[X] >> 4, sy = cal_step_size[Y] >> 4;
	 /*
  Serial.print("sx = ");
  Serial.println(sx);
  Serial.print("sy = ");
  Serial.println(sy);  
  */
	unsigned int x = val[X] - minVal[X];
	unsigned int y = val[Y] - minVal[Y];

	unsigned int xi0 = x / cal_step_size[X],
		yi0 = y / cal_step_size[Y];
	unsigned int xi1, yi1;
	if (xi0 * cal_step_size[X] >= size[X]) {
		xi1 = xi0;
		xi0 = xi1 - 1;
	}
	else {
		xi1 = xi0 + 1;
	}
	if (yi0 * cal_step_size[Y] >= size[Y]) {
		yi1 = yi0;
		yi0 = yi1 - 1;
	}
	else {
		yi1 = yi0 + 1;
	}
	unsigned long x0 = xi0 * sx,
		y0 = yi0 * sy,
		x1 = xi1 * sx,
		y1 = yi1 * sy;

	/*
	SERIAL_ECHOPAIR("x0 = ", x0);
	SERIAL_ECHOPAIR(" x1 = ", x1);
	SERIAL_ECHOPAIR(" y0 = ", y0);
	SERIAL_ECHOPAIR(" y1 = ", y1);
	SERIAL_ECHOLN("");
	*/

	static int q00[2], q01[2], q10[2], q11[2];
	static long a0, a1, a2, a3, a4, afx, afy;

	//loads the values of the 4 calibration points around the point requested
	q00[X] = offsets[xi0][yi0][X];
	q01[X] = offsets[xi0][yi1][X];
	q10[X] = offsets[xi1][yi0][X];
	q11[X] = offsets[xi1][yi1][X];

	q00[Y] = offsets[xi0][yi0][Y];
	q01[Y] = offsets[xi0][yi1][Y];
	q10[Y] = offsets[xi1][yi0][Y];
	q11[Y] = offsets[xi1][yi1][Y];

	// Really ugly bilinear equation from wikipedia...
	// Not fast, but currently only used during planner phase
	x = x >> 4;
	y = y >> 4;
	a0 = (x1 - x) * (y1 - y);
	a1 = (x - x0) * (y1 - y);
	a2 = (x1 - x) * (y - y0);
	a3 = (x - x0) * (y - y0);
	a4 = (x1 - x0) * (y1 - y0);
	afx = (q00[X] * a0 + q10[X] * a1 + q01[X] * a2 + q11[X] * a3) / a4;
	afy = (q00[Y] * a0 + q10[Y] * a1 + q01[Y] * a2 + q11[Y] * a3) / a4;
	/*
	SERIAL_ECHOPAIR("a0 = ", (float)a0);
	SERIAL_ECHOPAIR(" a1 = ", (float)a1);
	SERIAL_ECHOPAIR(" a2 = ", (float)a2);
	SERIAL_ECHOPAIR(" a3 = ", (float)a3);
	SERIAL_ECHOPAIR(" a4 = ", (float)a4);
	SERIAL_ECHOPAIR(" afx = ", (float)afx);
	*/
	val[X] += afx;
	val[Y] += afy;
}

void Galvo::ApplySlopeOffsets(unsigned int * val) {
	unsigned int x = val[X] - minVal[X];
	unsigned int y = val[Y] - minVal[Y];
	unsigned char index[2], index0[2];
	unsigned int origin[2];
	int diff[2];
	index[X] = select_index(x, X);
	index[Y] = select_index(y, Y);
	unsigned char t = 0,
				tx = index[X] < steps >> 1,
				ty = index[Y] < steps >> 1;
	MultiSU16X8toL16(origin[X], cal_step_size[X], index[X] + tx);//origin[X] = index[X] * cal_step_size[X];
	MultiSU16X8toL16(origin[Y], cal_step_size[Y], index[Y] + ty);//origin[Y] = index[Y] * cal_step_size[Y];
	index0[X] = index[X] + tx;
	index0[Y] = index[Y] + ty;
	/*
	if (index[X] < steps >> 1) {
		if (index[Y] < steps >> 1) {
			origin[X] = (index[X] + 1) * cal_step_size[X];//select_high_index_value(x, X);
			origin[Y] = (index[Y] + 1) + cal_step_size[Y];//select_high_index_value(y, Y);
			index0[X] = index[X] + 1;
			index0[Y] = index[Y] + 1;
		}
		else {
			origin[X] = (index[X] + 1) * cal_step_size[X];
			origin[Y] = index[Y] * cal_step_size[X];
			index0[X] = index[X] + 1;
			index0[Y] = index[Y];
		}
	}
	else {
		if (index[Y] < steps >> 1) {
			origin[X] = index[X] * cal_step_size[X];
			origin[Y] = (index[Y] + 1) * cal_step_size[X];
			index0[X] = index[X];
			index0[Y] = index[Y] + 1;
		}
		else {
			MultiSU16X8toL16(origin[X], cal_step_size[X], index[X]);//origin[X] = index[X] * cal_step_size[X];
			MultiSU16X8toL16(origin[Y], cal_step_size[Y], index[Y]);//origin[Y] = index[Y] * cal_step_size[Y];
			index0[X] = index[X];
			index0[Y] = index[Y];
		}
	}
	*/
	/*
	Serial.print("\nindex = ");
	printPair(index);
	Serial.print("\norigin = ");
	printPair(origin);
	Serial.print("\ndiff = ");
	printPair(diff);
	Serial.print("\nSlope X = ");
	printPair(slopes[index[X]][index[Y]][t][X]);
	Serial.print(" ");
	Serial.print((float)slopes[index[X]][index[Y]][t][X][X] / (float)FIXED_DIV);
	Serial.print(", ");
	Serial.print((float)slopes[index[X]][index[Y]][t][X][Y] / (float)FIXED_DIV);
	Serial.print("\nSlope Y = ");
	printPair(slopes[index[X]][index[Y]][t][Y]);
	Serial.println("");
	*/
	diff[X] = x - origin[X];
	diff[Y] = y - origin[Y];
	t = !(abs(diff[X] * cal_step_size[Y]) > abs(diff[Y] * cal_step_size[X]));
		/*
	if (abs(diff[X] * cal_step_size[Y]) > abs(diff[Y] * cal_step_size[X])) {
		t = 0;
	}
	else {
		t = 1;
	}
	*/
	long tmpDiff[2];
	MultiS16X16to32(tmpDiff[X], diff[X], slopes[index[X]][index[Y]][t][X][X]);
	MultiS16X16to32(tmpDiff[Y], diff[Y], slopes[index[X]][index[Y]][t][X][Y]);
	val[X] += (int)((tmpDiff[X] + tmpDiff[Y]) >> FIXED_BITS) + offsets[index0[X]][index0[Y]][X];
	MultiS16X16to32(tmpDiff[X], diff[X], slopes[index[X]][index[Y]][t][Y][X]);
	MultiS16X16to32(tmpDiff[Y], diff[Y], slopes[index[X]][index[Y]][t][Y][Y]);
	val[Y] += (int)((tmpDiff[X] + tmpDiff[Y]) >> FIXED_BITS) + offsets[index0[X]][index0[Y]][Y];
}

//Calculates an absolute galvo position based on a float input
void Galvo::CalcGalvoPosition(unsigned int * val, float x, float y) {
	//sane values
	x = max(min(x, mm_size[X]), 0.0);
	y = max(min(y, mm_size[Y]), 0.0);
	val[X] = minVal[X] + size[X] * x / mm_size[X];
	val[Y] = minVal[Y] + size[Y] * y / mm_size[Y];
	//unsigned long elapsed = micros();
	ApplyOffsets(val);
	//elapsed = micros() - elapsed;
	//SERIAL_ECHOPAIR("apply_offset elapsed time: ", elapsed);
	//SERIAL_ECHOLN("us");
}

/*
void Galvo::printCalTable() {
	for (int j = points - 1; j >= 0; j--) {
		for (int i = 0; i < points; i++) {
			printPair(offsets[i][j]);
		}
		Serial.println("");
	}
}



void Galvo::printValues() {
	Serial.print("min_step_size = ");
	printPair(min_step_size);
	Serial.print("\nmax_steps_per_unit = ");
	printPair(max_steps_per_unit);
	Serial.print("\nsteps = ");
	Serial.println(steps);
	Serial.print("points = ");
	Serial.println(points);
	Serial.print("scale = ");
	printPair(scale);
	Serial.print("\ncenter = ");
	printPair(center);
	Serial.print("\nmin = ");
	printPair(minVal);
	Serial.print("\nmax = ");
	printPair(maxVal);
	Serial.print("\nsize = ");
	printPair(size);
	Serial.print("\nmm_size = ");
	printPair(mm_size);
	Serial.print("\nsteps_per_mm = ");
	printPair(steps_per_mm);
	Serial.print("\ncal_step_size = ");
	printPair(cal_step_size);
	Serial.print("\ncal_step_size_mm = ");
	printPair(cal_step_size_mm);
	Serial.print("\nz_size = ");
	printPair(z_size);
	Serial.print("\nt0_max = ");
	printPair(t0_max);
	Serial.print("\nt_max = ");
	printPair(t_max);
	Serial.print("\ne = ");
	Serial.println(e);
	}
	*/
