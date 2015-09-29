#include "galvo.h"
#include "arduino.h"

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

//Calculates an absolute galvo position based on a float input
void Galvo::CalcGalvoPosition(unsigned int * val, float x, float y) {
	//sane values
	x = max(min(x, mm_size[X]), mm_size[X]);
	y = max(min(y, mm_size[Y]), mm_size[Y]);
	unsigned int x_tmp = minVal[X] + size[X] * x / mm_size[X];
	unsigned int y_tmp = minVal[Y] + size[Y] * y / mm_size[Y];
	val[X] = x_tmp;
	val[Y] = y_tmp;
	//unsigned long elapsed = micros();
	ApplyOffsets(val);
	//elapsed = micros() - elapsed;
	//SERIAL_ECHOPAIR("apply_offset elapsed time: ", elapsed);
	//SERIAL_ECHOLN("us");
}

void Galvo::printCalTable() {
	for (int j = 0; j < points; j++) {
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
	Serial.print("\nstep_size = ");
	printPair(step_size);
	Serial.print("\nz_size = ");
	printPair(z_size);
	Serial.print("\nt0_max = ");
	printPair(t0_max);
	Serial.print("\nt_max = ");
	printPair(t_max);
	Serial.print("\ne = ");
	Serial.println(e);
	}
// disabled cpp file as galvo has been moved to a class
#if 0
// Galvo related global variables
const float min_step_size = X_MAX_LENGTH / (float)0xFFFF;
const float max_steps_per_unit = 1.0 / min_step_size;
unsigned const int steps = CAL_GRID_SIZE;
unsigned const int points = steps + 1;
const unsigned int center = 0xFFFF / 2;

//Structure to hold an x/y offset per calibration point
offset offsets[points][points];

// The scale of the full grid with respect to full DAC space (0x0000 to 0xFFFF)
// This is used to initially populate x_min and x_max
float g_scale[2] = {
	GALVO_X_SCALE,
	GALVO_Y_SCALE
};

//The shift value for DAC space.  This will induce a tilt on an axis and used to compensate for a non vertical center point
unsigned int t_shift[2] = { 0, 0 };

//The center of the full grid with respect to full DAC space (0x0000 to 0xFFFF)
unsigned int g_center[2] = {
	center + t_shift[X_AXIS],
	center + t_shift[Y_AXIS]
};

//The calculated minimum DAC value
unsigned int g_min[2] = {
	(int)((float)g_center[X_AXIS] * (1 - g_scale[X_AXIS])),
	(int)((float)g_center[Y_AXIS] * (1 - g_scale[Y_AXIS]))
};

//The calculated maximum DAC value.  -1 is needed in order to make it an even value (hacky... may need some better math)
unsigned int g_max[2] = {
	(int)((float)g_center[X_AXIS] * (1 + g_scale[X_AXIS])) - 1,
	(int)((float)g_center[Y_AXIS] * (1 + g_scale[Y_AXIS])) - 1
};

//The calculated total size of the printable space in DAC units
unsigned int g_size[2] = {
	g_max[X_AXIS] - g_min[X_AXIS],
	g_max[Y_AXIS] - g_min[Y_AXIS]
};

//The distance between each calibration point
unsigned int cal_step_size[2] = {
	g_size[X_AXIS] / steps,
	g_size[Y_AXIS] / steps
};

// The calculated steps per mm (unrelated to printing steps per unit)
// This is used to calculate z_size and e
float steps_per_mm[2] = {
						(float)g_size[X_AXIS] / (float)X_MAX_POS,
						(float)g_size[Y_AXIS] / (float)Y_MAX_POS
						};

//The distance between the print bed and last galvo mirror expressed in steps
unsigned int z_size[2] = {
	(unsigned int)(Z_MAX_POS * steps_per_mm[X_AXIS]),
	(unsigned int)(Z_MAX_POS * steps_per_mm[Y_AXIS])
};

//The max theta (per axis) of the galvo as calculated by the size of the print area
const float t0_max = atan((X_MAX_POS / 2) / Z_MAX_POS);

//The max theta expressed in DAC steps
float t_max[2] = {
	(float)(g_size[X_AXIS] / 2.0) / t0_max,
	(float)(g_size[Y_AXIS] / 2.0) / t0_max
};

// The distance between the mirrors in DAC steps.  Using X axis as this distance 
// is the imaginary extension of the Z axis between the second mirror and the first
// which means it would be part of the X/Z coordinate system
const unsigned int e = E_DISTANCE * steps_per_mm[X_AXIS];

//Global galvo position. 
unsigned int g_position[2] = { g_min[X_AXIS], g_min[Y_AXIS] };

//Computes the calibration offset table
void compute_calibration_offsets() {
	float x_pos = (float)g_min[X_AXIS], y_pos = (float)g_min[Y_AXIS];
	float x_tmp, y_tmp;
	SERIAL_ECHO_START;
	SERIAL_ECHOLN("Calibration Table");
	/*
	SERIAL_ECHOPAIR("g_center = ", (unsigned long)g_center[Y_AXIS]);
	SERIAL_ECHOLN("");
	SERIAL_ECHOPAIR("g_min = ", (unsigned long)g_min[Y_AXIS]);
	SERIAL_ECHOLN("");
	SERIAL_ECHOPAIR("g_max = ", (unsigned long)g_max[Y_AXIS]);
	SERIAL_ECHOLN("");
	SERIAL_ECHOPAIR("g_size = ", (unsigned long)g_size[Y_AXIS]);
	SERIAL_ECHOLN("");
	SERIAL_ECHOPAIR("z_size = ", (unsigned long)z_size[Y_AXIS]);
	SERIAL_ECHOLN("");
	SERIAL_ECHOPAIR("t_max = ", (unsigned long)t_max[Y_AXIS]);
	SERIAL_ECHOLN("");
	*/
	for (int j = 0; j < points; j++) {
		for (int i = 0; i < points; i++) {
			x_tmp = x_pos - (float)g_center[X_AXIS];
			y_tmp = y_pos - (float)g_center[Y_AXIS];
			y_tmp = atan2(y_tmp, z_size[Y_AXIS]) * t_max[Y_AXIS]; // compute y first as its designated the independent (second galvo) axis
			x_tmp = atan2(x_tmp, sqrt((y_tmp * y_tmp) + ((float)z_size[X_AXIS] * (float)z_size[X_AXIS])) + (float)e) * t_max[X_AXIS];  // compute x next as its dependent on the Y movement
			x_tmp += (float)g_center[X_AXIS] - x_pos;
			y_tmp += (float)g_center[Y_AXIS] - y_pos;
			offsets[i][j].x = (int)x_tmp;
			offsets[i][j].y = (int)y_tmp;
			SERIAL_ECHO("(");
			SERIAL_ECHO(offsets[i][j].x);
			SERIAL_ECHO(", ");
			SERIAL_ECHO(offsets[i][j].y);
			SERIAL_ECHO(")");
			x_pos += cal_step_size[X_AXIS];
			if (x_pos > g_max[X_AXIS]) x_pos = g_min[X_AXIS];
		}
		y_pos += cal_step_size[Y_AXIS];
		if (y_pos > g_max[Y_AXIS]) y_pos = g_min[Y_AXIS];
		SERIAL_ECHOLN("");
	}
}

//Applies the offset table to a set of coordinates
void apply_offset(coord * val) {
	// shifting all coordinate values by 4 (divide by 16) in order to stop overflowing 32 bit longs
	// This reduces bilinear interpolation resolution, but is necessary in order to not use floats.  
	// This is not perfect and will still overflow if low values are chosen for number of grid points.

	unsigned int sx = cal_step_size[X_AXIS] >> 4, sy = cal_step_size[Y_AXIS] >> 4;
	unsigned int x = val->x - g_min[X_AXIS];
	unsigned int y = val->y - g_min[Y_AXIS];
	unsigned int xi0 = x / cal_step_size[X_AXIS],
		yi0 = y / cal_step_size[Y_AXIS];
	unsigned int xi1, yi1;
	if (xi0 * cal_step_size[X_AXIS] >= g_size[X_AXIS]) {
		xi1 = xi0;
		xi0 = xi1 - 1;
	}
	else {
		xi1 = xi0 + 1;
	}
	if (yi0 * cal_step_size[Y_AXIS] >= g_size[Y_AXIS]) {
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
	static offset q00, q01, q10, q11;
	static long a0, a1, a2, a3, a4, afx, afy;

	//loads the values of the 4 calibration points around the point requested
	q00.x = offsets[xi0][yi0].x;
	q01.x = offsets[xi0][yi1].x;
	q10.x = offsets[xi1][yi0].x;
	q11.x = offsets[xi1][yi1].x;

	q00.y = offsets[xi0][yi0].y;
	q01.y = offsets[xi0][yi1].y;
	q10.y = offsets[xi1][yi0].y;
	q11.y = offsets[xi1][yi1].y;

	// Really ugly bilinear equation from wikipedia...
	// Not fast, but currently only used during planner phase
	x = x >> 4;
	y = y >> 4;
	a0 = (x1 - x) * (y1 - y);
	a1 = (x - x0) * (y1 - y);
	a2 = (x1 - x) * (y - y0);
	a3 = (x - x0) * (y - y0);
	a4 = (x1 - x0) * (y1 - y0);
	afx = (q00.x*a0 + q10.x*a1 + q01.x*a2 + q11.x*a3) / a4;
	afy = (q00.y*a0 + q10.y*a1 + q01.y*a2 + q11.y*a3) / a4;
	/*
	SERIAL_ECHOPAIR("a0 = ", (float)a0);
	SERIAL_ECHOPAIR(" a1 = ", (float)a1);
	SERIAL_ECHOPAIR(" a2 = ", (float)a2);
	SERIAL_ECHOPAIR(" a3 = ", (float)a3);
	SERIAL_ECHOPAIR(" a4 = ", (float)a4);
	SERIAL_ECHOPAIR(" afx = ", (float)afx);
	*/
	val->x += afx;
	val->y += afy;
}

// calculates the absolute galvo position based on all offsets and calibration points
void abs_galvo_position(coord * val, float x, float y) {
	unsigned int x_tmp = g_min[X_AXIS] + g_size[X_AXIS] * x / max_pos[X_AXIS];
	unsigned int y_tmp = g_min[Y_AXIS] + g_size[Y_AXIS] * y / max_pos[Y_AXIS];
	val->x = x_tmp;
	val->y = y_tmp;
	//unsigned long elapsed = micros();
	apply_offset(val);
	//elapsed = micros() - elapsed;
	//SERIAL_ECHOPAIR("apply_offset elapsed time: ", elapsed);
	//SERIAL_ECHOLN("us");
}

float step_size[2] = { X_MAX_LENGTH / CAL_GRID_SIZE, Y_MAX_LENGTH / CAL_GRID_SIZE };
float get_x(int i) { return X_MIN_POS + step_size[X_AXIS] * i; }
float get_y(int i) { return Y_MIN_POS + step_size[Y_AXIS] * i; }

int select_x_index(float x) {
	int i = 1;
	while (x > get_x(i) && i < steps) i++;
	return i - 1;
}

int select_y_index(float y) {
	int i = 1;
	while (y > get_y(i) && i < steps) i++;
	return i - 1;
}
#endif