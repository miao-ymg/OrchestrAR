#include "../include/Instrument.h"
#include <math.h>

#include <random>

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
	
	std::array<float, 3> rgb = paintObject(pitch);
	float r = rgb[0];
	float g = rgb[1];
	float b = rgb[2];
	glColor4f(r, g, b, 1);

	switch (role) {
		case BASS:
			drawSpeaker(0.5, 0.5, 1, r, g, b);
			break;
		case BEAT:
            drawDrum(0.5, 0.5, r, g, b);
			break;
		case KEYS:
			drawPianoKeys(0.9, 0.25, 0.15, r, g, b);
			break;
		case MELODY:
			drawMusicalNote();
			glTranslatef(0.35, 0, 0.2);
			drawMusicalNote();
			glTranslatef(-0.7, 0, 0.15);
			drawMusicalNote();
			break;
		case VOCAL:
            drawMicrophone(r, g, b);
			break;
		default:
            throw std::invalid_argument("ERROR: Sample role shouldn't exist!");
	}
	
}

void Instrument::freeChunk() {
    Mix_FreeChunk(sample);
}
