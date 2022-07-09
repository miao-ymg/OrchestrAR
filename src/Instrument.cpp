#include "../include/Instrument.h"


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


void Instrument::toggleSound(vector<int>& identifiers) {
    if (find(identifiers.begin(), identifiers.end(), id) != identifiers.end()) {
		Mix_Volume(channel, volume);
		// std::cout << "Vol: " << volume << std::endl;
	} else {
		--timeToLive;
		if (timeToLive == 0)
			Mix_Volume(channel, 0);
	}

	// --- Repeat the sound file if neccesary ---
	if (Mix_Playing(channel) != 1)
        Mix_PlayChannel(channel, sample, 0);
}


void Instrument::freeChunk() {
    Mix_FreeChunk(sample);
}
