#include "../include/Instrument.h"
#include <math.h>

#include <random>
//#include "DrawPrimitives.h"

/* PI */
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

using namespace std;

const int MUSIC_STOP_DELAY = 16;

Instrument::Instrument(const int id, Role role, Pitch pitch, const char* soundPath) 
	: 
	soundPath(std::move(soundPath)), id(id), timeToLive(MUSIC_STOP_DELAY), role(role), pitch(pitch) {};


int Instrument::getID() const { return id; }


Role Instrument::getRole() const { return role; }


Pitch Instrument::getPitch() const { return pitch; }


int Instrument::getVolume() const { return volume; }


std::array<float, 16> Instrument::getPoseMatrix() {
	std::array<float, 16> ret;
	for (size_t i = 0; i < 16; i++)
		ret[i] = poseMatrix[i];
	return ret;
}


void Instrument::setPoseMatrix(float matrix[16]){
	memcpy(poseMatrix, matrix, sizeof(float) * 16);
}


void Instrument::setVolume(double volume) {
	this->timeToLive = MUSIC_STOP_DELAY;
	this->volume = static_cast<int>(volume * 128);
}


void Instrument::loadSound() {
    sample = Mix_LoadWAV(soundPath);

	if (sample == NULL)
		fprintf(stderr, "Unable to load WAV file: %s\n", Mix_GetError());

	channel = Mix_PlayChannel(-1, sample, 0);
	if (channel == -1)
		fprintf(stderr, "Unable to play WAV file: %s\n", Mix_GetError());

	Mix_Volume(channel, 0);
}

void Instrument::loadRandomSound(Role role) {
	loadSound();
}

void Instrument::toggleSound(vector<int>& identifiers) {
    if (find(identifiers.begin(), identifiers.end(), id) != identifiers.end()) {
		Mix_Volume(channel, volume);
	} else {
		--timeToLive;
		if (timeToLive == 0)
			Mix_Volume(channel, 0);
	}

	// --- Repeat the sound file if neccesary ---
	if (Mix_Playing(channel) != 1)
        Mix_PlayChannel(channel, sample, 0);
}


void Instrument::drawObject(){
	glTranslatef(0, 0, 0.1);
	// Make alpha value depending on volume (between 0 and 1)
	// float alpha = (float) instr.getVolume() / 128;
	float alpha = 1;
	/*
	float r, g, b;
	switch (role) {
		case BASS:
			r = 0.5; g = 0.5; b = 0.5;
			break;
		case BEAT:
			r = 0; g = 1; b = 0;
			break;
		case KEYS:
			r = 0; g = 0.5; b = 0.5;
			break;
		case MELODY:
			r = 0; g = 0; b = 1;
			break;
		case VOCAL:
			r = 0; g = 1; b = 0;
			break;
		default:
			r = 0; g = 1; b = 1;
			break;
	}
	*/
	std::array<float, 3> rgb = setColor();
	float r = rgb[0];
	float g = rgb[1];
	float b = rgb[2];
	glColor4f(r, g, b, alpha);

	float drumWidth = 0.5;
	float drumHeight = 0.5;

	switch (role) {
		case BASS:
			drawSpeakers(0.5, 0.5, 1, r, g, b);
			break;
		case BEAT:
			for (size_t i = 0; i < 6; i++) {
				glTranslatef(drumWidth * 2 * M_PI * sin(i) / 6, drumWidth * 2 * M_PI * cos(i) / 6, 0);
				drawCylinder(0.03, 0.7, 1, 1, 1);
				glTranslatef(drumWidth * -2 * M_PI * sin(i) / 6, drumWidth * -2 * M_PI * cos(i) / 6, 0);
			}
			// Drum shell
			glTranslatef(0, 0, 0.2);
			drawCylinder(drumWidth, drumHeight, r, g, b);
			// Drum skin
			glTranslatef(0, 0, drumHeight - 0.0);
			drawCylinder(drumWidth + 0.05, 0.05, 1, 1, 1);
			break;
		case KEYS:
			drawPianoKeys(1, 0.2, 0.15);
			break;
		case MELODY:
			drawMusicNote();
			glTranslatef(0.35, 0, 0.2);
			drawMusicNote();
			glTranslatef(-0.7, 0, 0.15);
			drawMusicNote();
			break;
		case VOCAL:
			glRotatef(15, 1, 0, 0);
			drawCylinder(0.15, 1.5, r, g, b);
			glTranslatef(0, 0, 1.4);
			drawCylinder(0.3, 0.04, r, g, b);
			glColor4f(0.2, 0.2, 0.2, alpha);
			drawSphere(0.3, 10, 10);
			glTranslatef(0, 0, -2);
			break;
		default:
			break;
	}
	
}

void Instrument::freeChunk() {
    Mix_FreeChunk(sample);
}

std::array<float, 3> Instrument::setColor() {
	if (pitch >= 24)
		return {1.0, 25.0/255, 179.0/255};

	float r = 0.5 * (cos(2.0 * M_PI * pitch / 24) + 1);
	float g = 0.5 * (cos(2.0 * M_PI * (pitch - 8) / 24) + 1);
	float b = 0.5 * (cos(2.0 * M_PI * (pitch + 8) / 24) + 1);

	return {r, g, b};
}