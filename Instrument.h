#pragma once


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
    void playSound(vector<int>& identifiers);

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