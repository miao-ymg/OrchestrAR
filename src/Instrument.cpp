
#include <vector>
#include <algorithm>

#include <SDL.h>
#include <SDL_mixer.h>

#include "../include/Instrument.h"


using namespace std;

Instrument::Instrument(const int id, const char* soundPath) 
	: 
	soundPath(std::move(soundPath)), id(id) {};


int Instrument::getID() const {
	return id;
}


void Instrument::setVolume(double volume) {
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
		std::cout << "Vol: " << volume << std::endl;
	} else {
		Mix_Volume(channel, 0);
	}

	// --- Repeat the sound file if neccesary ---
	if (Mix_Playing(channel) != 1)
        Mix_PlayChannel(channel, sample, 0);
}


void Instrument::freeChunk() {
    Mix_FreeChunk(sample);
}
