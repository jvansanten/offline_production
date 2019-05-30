
#include "HepRandomEngineFromI3RandomService.h"

HepRandomEngineFromI3RandomService::HepRandomEngineFromI3RandomService(I3RandomServicePtr rng) : randomService_(rng)
{}

HepRandomEngineFromI3RandomService::~HepRandomEngineFromI3RandomService() {}

double HepRandomEngineFromI3RandomService::flat()
{
    return randomService_->Uniform();
}

void HepRandomEngineFromI3RandomService::flatArray(const int size, double *vect)
{
    for (int i=0; i<size; i++)
        vect[i] = randomService_->Uniform();
}

