#pragma once
#ifndef GALVO_H
#define GALVO_H
//#include "PDQ_FastPin.h"
#include <math.h>
#include "arduino.h"

#define TRIANGLE_MESH
#ifdef TRIANGLE_MESH
#define FIXED_BITS 14
#define FIXED_DIV 16384
#endif

#define CAL_GRID_SIZE 2
class Galvo {
public:
	enum {
		X = 0,
		Y = 1,
		Z = 3,
		E = 4,
		S0 = 0,
		S1 = 1,
		DAC_MAX = 0xFFFF
	};
private:
	float min_step_size[2];
	float max_steps_per_unit[2];
	const static unsigned int steps = CAL_GRID_SIZE;
	const static unsigned int points = CAL_GRID_SIZE + 1;

	//The scale of the full grid with respect to full DAC space (0x0000 to 0xFFFF)
	float scale[2];

	// Max value of the DAC.  Could concievably need to to be modified so it exists as a variable
	unsigned int maxDAC = 0xFFFF;

	//The shift value for DAC space.  This will induce a tilt on an axis and used to compensate for a non vertical center point
	unsigned int shift[2] = { 0, 0 };

	//The center of the full grid with respect to full DAC space (0x0000 to 0xFFFF)
	unsigned int center[2] = {
		maxDAC/2 + shift[X],
		maxDAC/2 + shift[Y]
	};

	//The calculated minimum DAC value
	unsigned int minVal[2];

	//The calculated maximum DAC value.  -1 is needed in order to make it an even value (hacky... may need some better math)
	unsigned int maxVal[2];

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
	float cal_step_size_mm[2];

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
	//array to hold table of slopes for mesh/slope calculations
	int slopes[steps][steps][2][2][2] = {};
#endif

	Galvo();

	Galvo(float x_length, float y_length, float z_length, float e_length, float x_scale, float y_scale);

	// Initializes all variables
	// This should be called after you vary any initialized parameter
	void Initialize();

	//Computes the calibration offset table
	void CalcCalibrationTable();

	//Zeroes the calibration offset table
	void ResetCalibrationTable();

#ifdef TRIANGLE_MESH
	//Calculates the slope table in order to use triangle mesh
	void CalcSlopeTable();
	//Zeroes the slope table
	void ResetSlopeTable();
#endif

	//Applies the offset table to a set of coordinates
	void ApplyOffsets(unsigned int * val);

	//Applies the offset table to a set of coordinates using the slope table
	void ApplySlopeOffsets(unsigned int * val);

	// calculates the absolute galvo position based on all offsets and calibration points
	void CalcGalvoPosition(unsigned int * val, float x, float y);

	// gets the grid index 
	float get_x_mm(int i) { return cal_step_size_mm[X] * i; }
	float get_y_mm(int i) { return cal_step_size_mm[Y] * i; }

	int select_x_index(float x) {
		int i = 1;
		while (x > get_x_mm(i) && i < steps) i++;
		return i - 1;
	}

	int select_y_index(float y) {
		int i = 1;
		while (y > get_y_mm(i) && i < steps) i++;
		return i - 1;
	}

	inline int select_index(unsigned int val, unsigned char axis) {
		int i = 1;
		unsigned int val0 = 0;
		while (val > val0 && i < steps) {
			i++;
			val0 += cal_step_size[axis];
		}
		return i - 1;
	}

	inline int select_low_index_value(unsigned int val, unsigned char axis) {
		int i = 1;
		unsigned int val0 = cal_step_size[axis];
		while (val > val0 && i < steps) {
			i++;
			val0 += cal_step_size[axis];
		}
		return val0 - cal_step_size[axis];
	}

	inline int select_high_index_value(unsigned int val, unsigned char axis) {
		int i = 1;
		unsigned int val0 = cal_step_size[axis];
		while (val > val0 && i < steps) {
			i++;
			val0 += cal_step_size[axis];
		}
		return val0;
	}

	void calcMinStepSize() {
		min_step_size[X] = mm_size[X] / (float)maxDAC;
		min_step_size[Y] = mm_size[Y] / (float)maxDAC;
	}

	void calcMaxStepsPerUnit() {
		max_steps_per_unit[X] = 1.0 / min_step_size[X];
		max_steps_per_unit[Y] = 1.0 / min_step_size[Y];
	}

	void calcMin() {
		minVal[X] = (int)((float)center[X] * (1.0 - scale[X]));
		minVal[Y] = (int)((float)center[Y] * (1.0 - scale[Y]));
	};

	void calcMax() {
		maxVal[X] = (int)((float)center[X] * (1.0 + scale[X])) - 1;
		maxVal[Y] = (int)((float)center[Y] * (1.0 + scale[Y])) - 1;
	};

	void calcSize() {
		size[X] = maxVal[X] - minVal[X];
		size[Y] = maxVal[Y] - minVal[Y];
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

	void calcStepSize() {
		cal_step_size_mm[X] = mm_size[X] / (float)steps;
		cal_step_size_mm[Y] = mm_size[Y] / (float)steps;
	}

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

	unsigned int getMin(unsigned char axis) {
		return minVal[axis];
	}

	unsigned int getMax(unsigned char axis) {
		return maxVal[axis];
	}

	float getStepSize(unsigned char axis) {
		return cal_step_size_mm[axis];
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

	unsigned int getMaxStepsPerUnit(unsigned char axis) {
		return max_steps_per_unit[axis];
	}
#ifdef GALVO_DEBUG
	template <typename T> static void printPair(T * val) {
		Serial.print("(");
		Serial.print(val[X]);
		Serial.print(", ");
		Serial.print(val[Y]);
		Serial.print(")");
	}
	
	void printValues();

	void printCalTable();

	void printSlopeTable();
#endif
};

#endif