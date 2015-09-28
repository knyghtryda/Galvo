#pragma once
#ifndef GALVO_H
#define GALVO_H
//#include "PDQ_FastPin.h"
#include <math.h>

class Galvo {
private:
	enum {
		X = 0,
		Y = 1,
		Z = 3,
		E = 4,
		CAL_GRID_SIZE = 8,
		DAC_MAX = 0xFFFF
	};
	float min_step_size;
	float max_steps_per_unit;
	const static unsigned int steps = CAL_GRID_SIZE;
	const static unsigned int points = CAL_GRID_SIZE + 1;

	//The scale of the full grid with respect to full DAC space (0x0000 to 0xFFFF)
	float scale[2];

	// number of bits in the DAC.  Could concievably need to to be modified so it exists as a variable
	unsigned int dacBits = DAC_MAX;

	//The shift value for DAC space.  This will induce a tilt on an axis and used to compensate for a non vertical center point
	unsigned int shift[2] = { 0, 0 };

	//The center of the full grid with respect to full DAC space (0x0000 to 0xFFFF)
	unsigned int center[2] = {
		dacBits/2 + shift[X],
		dacBits/2 + shift[Y]
	};

	//The calculated minimum DAC value
	unsigned int min[2];

	//The calculated maximum DAC value.  -1 is needed in order to make it an even value (hacky... may need some better math)
	unsigned int max[2];

	//The calculated total size of the printable space in DAC units
	unsigned int size[2];
	
	//The size of the target in mm
	float mm_size[4];

	// The calculated steps per mm (unrelated to printing steps per unit)
	// This is used to calculate z_size and e
	float steps_per_mm[2];

	//The distance between each calibration point
	unsigned int cal_step_size[2];

	//The distance between each calibration point in mm
	float step_size[2] = { mm_size[X] / steps, mm_size[Y] / steps };

	//The distance between the print bed and last galvo mirror expressed in steps
	unsigned int z_size[2];

	//The max theta (per axis) of the galvo as calculated by the size of the print area
	float t0_max[2];

	//The max theta expressed in DAC steps
	float t_max[2];

	//The distance between the mirrors in DAC steps
	unsigned int e;

	//Global galvo position. 
	unsigned int position[2];

public:
	//array to hold an x/y offset per calibration point
	int offsets[points][points][2] = {};

#ifdef TRIANGLE_MESH
	int slopes[steps][steps][3] = {};
#endif

	Galvo();

	Galvo(float x_length, float y_length, float z_length, float e_length, float x_scale, float y_scale);

	// Initializes all variables
	// This should be called after you vary any initialized parameter
	void Initialize();

	//Computes the calibration offset table
	void CalcCalibrationTable();

#ifdef TRIANGLE_MESH
	//Calculates the slope table in order to use triangle mesh
	void CalcSlopeTable();
#endif

	//Applies the offset table to a set of coordinates
	void ApplyOffsets(unsigned int * val);

	// calculates the absolute galvo position based on all offsets and calibration points
	void CalcGalvoPosition(unsigned int * val, float x, float y);

	// gets the grid index 
	float get_x(int i) { return min[X] + step_size[X] * i; }
	float get_y(int i) { return min[Y] + step_size[Y] * i; }

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

	void calcMin() {
		min[X] = (int)((float)center[X] * (1 - scale[X]));
		min[Y] = (int)((float)center[Y] * (1 - scale[Y]));
	};

	void calcMax() {
		max[X] = (int)((float)center[X] * (1 + scale[X])) - 1;
		max[Y] = (int)((float)center[Y] * (1 + scale[Y])) - 1;
	};

	void calcSize() {
		size[X] = max[X] - min[X];
		size[Y] = max[Y] - min[Y];
	};

	void setMMSize(float x, float y, float z, float e) {
		mm_size[X] = x;
		mm_size[Y] = y;
		mm_size[Z] = z;
		mm_size[E] = e;
	};

	void setScale(float x, float y) {
		scale[X] = x;
		scale[Y] = y;
	};
	
	void calcZSize() {
		z_size[X] = (unsigned int)(mm_size[Z] * steps_per_mm[X]);
		z_size[Y] = (unsigned int)(mm_size[Z] * steps_per_mm[Y]);
	};

	void calcCalStepSize() {
		cal_step_size[X] = size[X] / steps;
		cal_step_size[Y] = size[Y] / steps;
	};

	void calcT0Max() {
		t0_max[X] = atan((mm_size[X] / 2) / mm_size[Z]);
		t0_max[Y] = atan((mm_size[Y] / 2) / mm_size[Z]);
	};

	void calcTMax() {
		t_max[X] = (float)(size[X] / 2.0) / t0_max[X];
		t_max[Y] = (float)(size[Y] / 2.0) / t0_max[Y];
	};

	void calcStepsPerMM() {
		steps_per_mm[X] = (float)size[X] / (float)mm_size[X];
		steps_per_mm[Y] = (float)size[Y] / (float)mm_size[Y];
	};

	void setPosition(unsigned int x, unsigned int y) {
		position[X] = x;
		position[Y] = y;
	};

	void setOffset(unsigned int x, unsigned int y, int x_offset, int y_offset) {
		offsets[x][y][X] = x_offset;
		offsets[x][y][Y] = y_offset;
	}

	int getOffset(unsigned int x, unsigned int y, unsigned int axis) {
		return offsets[x][y][axis];
	}

	unsigned int getOffsetsSize() {
		return sizeof(offsets);
	}

	static const unsigned int getSteps() {
		return steps;
	}

	static const unsigned int getPoints() {
		return points;
	}

	unsigned int getMaxStepsPerUnit() {
		return max_steps_per_unit;
	}
};

#endif