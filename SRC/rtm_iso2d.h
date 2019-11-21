#ifndef RTM_ROUTINE_H
#define RTM_ROUTINE_H

#define VEC_LEN 1000

// -------- Temporary FILES...
static char crossCorrelationFile[] = "./TEMPORARY_FILES/crossCorrelation_%05dX%05d.bin";
static char selfSourceCorrelationFile[] = "./TEMPORARY_FILES/selfCrossCorrelation_%05dX%05d.bin";
// ----------------------------------

void rtm_routine(int simulationId);

#endif