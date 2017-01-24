#define ERR_INVALID_INPUT 1
#define ERR_MALLOC_FAIL 2
#define ERR_FILE_OPEN 3
#define ERR_PGPLOT 4


// newtypes
typedef struct {
    double x;
    double y;
    double z;
} vector;
typedef struct {
    vector r;
    vector rold;
    vector v;
    vector a;
    float m;
} particle;
typedef struct {
    float dt;
    float t;
    float G;
    int bounded;
    double lengthX;
    double lengthY;
    double lengthZ;
    int accretion;
    int nparticles;
    particle *p;
} environment;

void orbitEulerStep(environment *env, particle *p1, particle *p2);
void orbitRKStep(environment *env, particle *p1, particle *p2);
void orbitLeapfrogStep(environment *env);

float dvxdtGrav(float t, float x, float y);
float dvydtGrav(float t, float x, float y);

void doAccretion(environment *env);

void particleCombine(particle *p1, particle *p2);
void particleDeleteShift(environment *env, int n);
void particleSetR(particle *p, double x, double y, double z);
void particleSetV(particle *p, double x, double y, double z);
void particleSetVTangent(particle *p, double G, double M, double centerX, double centerY, double centerZ);
void initParticleZeroes(particle *p);
void initParticleRandom(particle *p);

// miscellaneous
void errorCase(const int errorCode);
