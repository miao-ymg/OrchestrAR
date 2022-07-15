#include "../include/Instrument.h"

#include <random>
//#include "DrawPrimitives.h"

using namespace std;

const int MUSIC_STOP_DELAY = 16;

Instrument::Instrument(const int id, Role role, const char* soundPath) 
	: 
	soundPath(std::move(soundPath)), id(id), timeToLive(MUSIC_STOP_DELAY), role(role) {};


int Instrument::getID() const { return id; }


Role Instrument::getRole() const { return role; }


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
	float r, g, b;
	switch (role) {
		case BASS:
			r = 1; g = 1; b = 1;
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
	glColor4f(r, g, b, alpha);

	float drumWidth = 0.5;
	float drumHeight = 0.5;

	switch (role) {
		case BASS:
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
			glTranslatef(0, 0, drumHeight - 0.04);
			drawCylinder(drumWidth + 0.05, 0.05, 1, 1, 1);
			break;
		case KEYS:
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
