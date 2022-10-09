#pragma once

#include <vector>
#include <SDL.h>
#include <SDL_mixer.h>


enum Role { BASS, BEAT, KEYS, MELODY, VOCAL };

enum Pitch {C, Em, G, Bm, D, FHm,
            A, CHm, E, GHm, B, DHm,
            FH, AHm, CH, Fm, Ab, Cm,
            Eb, Gm, Bb, Dm, F, Am,
            None };

class Instrument {
public:
    /**
     * Constructor
     */
    Instrument(const int id, Role role, Pitch pitch, const char* soundPath);
    
    /**
     * Getter for id
     */
    int getID() const;

    /**
     * Getter for role
     */
    Role getRole() const;

    /**
     * Getter for key
     */
    Pitch getPitch() const;

    /**
     * Getter for volume
     */
    int getVolume() const;

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
     * Load random sound file.
     */
    void loadRandomSound(Role role);

    /**
     * Play the sound in an endless loop
     */
    void toggleSound(std::vector<int>& identifiers);

    /**
     * Free the memory space occupied by the music chunk.
     */
    void freeChunk();


    void drawObject();

    /**
     * Convert a role to a string
     */
    static char* roleToString(Role r);


protected:
    const int id;
    const char* soundPath;
    const Role role;
    const Pitch pitch;

    Mix_Chunk* sample;
    float poseMatrix[16];
    int channel;
    int volume;
    int timeToLive;         // Remaining music playing time when marker disappears
};
