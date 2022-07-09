#pragma once

#include <vector>
#include <SDL.h>
#include <SDL_mixer.h>


class Instrument {
public:
    /**
     * Constructor
     */
    Instrument(const int id, const char* soundPath);
    
    /**
     * Getter for id
     */
    int getID() const;

    /**
     * Getter for poseMatrix
     */
    std::array<float, 16> getPoseMatrix();

    /**
     * Set pose informationi for this instruments marker
     */
    void setPoseMatrix(float matrix[16]);

    /**
     * Map the relative volume to absolute volume
     */
    void setVolume(double volume);

    /**
     * Load sound file.
     */
    void loadSound();

    /**
     * Play the sound in an endless loop
     */
    void toggleSound(std::vector<int>& identifiers);

    /**
     * Free the memory space occupied by the music chunk.
     */
    void freeChunk();


protected:
    const int id;
    const char* soundPath;

    Mix_Chunk* sample;
    float poseMatrix[16];
    int channel;
    int volume;
    int timeToLive;         // Remaining music playing time when marker disappears
};
