#pragma once

#include <vector>
#include <SDL.h>
#include <SDL_mixer.h>


using namespace std;

class Instrument {
public:
    /**
     * Constructor
     */
    Instrument(const int id, const char* soundPath);

    /**
     * Load sound file.
     */
    void loadSound();

    /**
     * Play the sound in an endless loop
     */
    void toggleSound(vector<int>& identifiers);

    /**
     * Free the memory space occupied by the music chunk.
     */
    void freeChunk();


protected:
    const int id;
    const char* soundPath;

    Mix_Chunk* sample;
    int channel;
};