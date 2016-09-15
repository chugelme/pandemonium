//
//  main.cpp
//  Pandemonium
//
//  Created by Cole Hugelmeyer on 9/8/16.
//  Copyright © 2016 DeadFish. All rights reserved.
//
#include <GLFW/glfw3.h>
#include <iostream>
#include <OpenGL/gl3.h>
#include <string.h>
#include <chrono>
#include <tgmath.h>
#include <math.h>
#include <cstdlib>
#define bodyparticles 100
#define auranum 100
#define BOUND 0
#define UNBOUND 1
#define DEAD 2
using namespace std;

class Monster1{
public:
    bool alive;
    float x;
    float y;
    float ang;
    float sang;
    float svel;
    float spe;
    float HP;
    float minorHP;
    float centerx;
    float centery;
    float centerx0;
    float centery0;
    float territoryRad;
    float sdir;
    float recharge[5];
    void increment(float bodx, float body, float vx, float vy, float dt);
    void move(float dx,float dy);
    void respawn();
};

void Monster1::increment(float bodx, float body, float vx, float vy, float dt){
    if(alive){
    if((bodx - centerx)*(bodx - centerx) + (body - centery)*(body - centery) > territoryRad*territoryRad){
        if((x - centerx)*(x - centerx) + (y - centery)*(y - centery) > 0.01*territoryRad){
            ang = 3.1415926 + atan2(y - centery,x - centerx);
            spe = 0.3;
            svel = 1;
            sdir = 1;
        }else{
            spe = 0;
            svel = 1;
        }
    }else{
        if((bodx - centerx)*(bodx - centerx) + (body - centery)*(body - centery) < territoryRad*territoryRad*0.1
           && (x - centerx)*(x - centerx) + (y - centery)*(y - centery) < territoryRad*territoryRad*0.5){
            ang = atan2(centery- (centery - body) - y, centerx- (centerx - bodx) - x) + 3.1415926;
            spe = 0.6;
            svel = 6*sdir;
        }else{
            ang = atan2(centery+ (centery - body)*0.8 - y, centerx+ (centerx - bodx)*0.8 - x);
            spe = 0.5;
            svel = 6*sdir;
        }
    }
    
    if (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) < dt*0.7){
        sdir = -sdir;
    }
    
    x += spe*cos(ang)*dt;
    y += spe*sin(ang)*dt;
    sang += svel*dt;
    
    for(int i = 0; i < 5; i++){
        if(recharge[i] < 1-(dt/1.7)){
            recharge[i] += dt/1.7;
        }else{
            recharge[i] = 1;
        }
        float dist = sqrt((body-y)*(body-y) + (bodx-x)*(bodx-x));
        if(abs( sin( (sang + (2* 3.1415926* ((float)i) )/5 - atan2(body + vy*dist/2 - y,bodx + vx*dist/2 - x) )/2 ) )< 0.03  && recharge[i] == 1){
            recharge[i] = 0;
        }
    }
        if(HP <= 0 && minorHP <=0){
            alive = false;
        }
    }
}

void Monster1::move(float dx, float dy){
    if(alive){
    x += dx;
    y += dy;
    centerx += dx;
    centery += dy;
    }
}

void Monster1::respawn(){
    alive = true;
    x = centerx0;
    y = centery0;
    centerx = centerx0;
    centery = centery0;
    HP = 1;
    minorHP = 0;
    for(int i = 0; i < 5; i++){
        recharge[i] = 1;
    }
    ang= 0;
    sang = 0;
    svel = 0;
    spe = 0;
    sdir = 1;
}

void IncrementParticles(float particles[], float parvel[], int N, float increment,
                        float xpos, float ypos, float dir, float dis, float ang, float vel){
    
    for(int i = 0; i < N; i++){
        particles[3*i + 2] -= 0.001;
        if(particles[3*i + 2] > 0){
            particles[3*i] += increment * parvel[2*i];
            particles[3*i + 1] += increment * parvel[2*i + 1];
        }else{
            float r = dis * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            float theta = 2 * 3.1415926 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5);
            particles[3*i] = r * cos(dir + theta) + xpos;
            particles[3*i + 1] = r * sin(dir + theta) + ypos;
            particles[3*i + 2] = 0.01 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            parvel[2*i] = vel * cos(dir - ang * theta);
            parvel[2*i + 1] = vel * sin(dir - ang * theta);
        }
    }
    
}

void MoveBody(float body[], int NB, float dx, float dy){
    for(int i = 0; i < NB; i++){
        body[3*i]+=dx;
        body[3*i + 1]+=dy;
    }
}

void MoveParticles(float particles[], int NP, float dx, float dy){
    for(int i = 0; i < NP; i++){
        particles[3*i]+=dx;
        particles[3*i + 1]+=dy;
    }
}

void IncrementBackground(float background[], float C, float S){
    float temp1;
    float temp2;
    float r;
    for(int i = 0; i < 9; i++){
        for(int j = 0; j < 9; j++){
            temp1 = background[4*(9*i + j) + 2];
            temp2 = background[4*(9*i + j) + 3];
            r = temp1*temp1 + temp2*temp2;
            background[4*(9*i + j) + 2] = (temp1*C - temp2*S) ;
            background[4*(9*i + j) + 3] = (temp2*C + temp1*S) ;
        }
    }
}


void MoveBackground(float background[], float dx, float dy, float density, float AspectRatio){
    for(int i = 0; i < 9*9; i++){
        background[4*i] += dx;
        background[4*i+1] += dy*AspectRatio;
    }
    
    // delete the following to make non-regenerating background
    
    
    if(background[9*4*4 + 4*4] > 1.7320508/density){
        for(int i = 0; i < 9*9; i++){
            background[4*i] -= 1.7320508/density;
        }

        for(int i = 8; i > 1; i--){
            for(int j = 0; j < 9; j++){
                for(int k = 2; k < 4;k++){
                    
                    background[4*(9*i + j) + k] = background[4*(9*(i-2) + j) + k];
                    
                }
            }
        }
        
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 9; j++){
                background[4*9*i + 4*j + 2] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
                background[4*9*i + 4*j + 3] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
            }
        }
        
    }
    if(background[9*4*4 + 4*4] < -1.7320508/density){
        for(int i = 0; i < 9*9; i++){
            background[4*i] += 1.7320508/density;
        }
        
        for(int i = 0; i < 7; i++){
            for(int j = 0; j < 9; j++){
                for(int k = 2; k < 4;k++){
                    
                    background[4*(9*i + j) + k] = background[4*(9*(i+2) + j) + k];
                    
                }
            }
        }
        
        for(int i = 7; i < 9; i++){
            for(int j = 0; j < 9; j++){
                background[4*9*i + 4*j + 2] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
                background[4*9*i + 4*j + 3] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
            }
        }
    }
    if(background[9*4*4 + 4*4 + 1] > 2*AspectRatio/density){
        for(int i = 0; i < 9*9; i++){
            background[4*i+1] -= 2*AspectRatio/density;
        }
        
        for(int i = 0; i < 9; i++){
            for(int j = 8; j > 1; j--){
                for(int k = 2; k < 4;k++){
                    
                    background[4*(9*i + j) + k] = background[4*(9*i + j-2) + k];
                    
                }
            }
        }
        
        for(int i = 0; i < 9; i++){
            for(int j = 0; j < 2; j++){
                background[4*9*i + 4*j + 2] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
                background[4*9*i + 4*j + 3] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
            }
        }
    }
    if(background[9*4*4 + 4*4 + 1] < - 2*AspectRatio/density){
        for(int i = 0; i < 9*9; i++){
            background[4*i+1] += 2*AspectRatio/density;
        }
        
        for(int i = 0; i < 9; i++){
            for(int j = 0; j < 7; j++){
                for(int k = 2; k < 4;k++){
                    
                    background[4*(9*i + j) + k] = background[4*(9*i + j+2) + k];
                    
                }
            }
        }
        
        for(int i = 0; i < 9; i++){
            for(int j = 7; j < 9; j++){
                background[4*9*i + 4*j + 2] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
                background[4*9*i + 4*j + 3] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
            }
        }

    }
}

void IncrementAura(float aura[], float auravel[], int auradata[], float aurarot, float dt, float bodx, float body, float mousex, float mousey, int regen){
    float temp1;
    float temp2;
    float C = cos(aurarot * dt);
    float S = sin(aurarot * dt);
    int count = regen;
    for(int i = 0; i < auranum; i++){
        aura[3*i] += dt* auravel[2*i];
        aura[3*i+1] += dt* auravel[2*i+1];
        if(auradata[i] == BOUND){
            temp1 = aura[3*i] - bodx;
            temp2 = (aura[3*i + 1]-body);
            aura[3*i] = C*temp1 - S*temp2 + bodx;
            aura[3*i+1] = (C*temp2 + S*temp1) + body;
            temp1 = auravel[2*i];
            temp2 = (auravel[2*i + 1]);
            auravel[2*i] = C*temp1 - S*temp2;
            auravel[2*i+1] = (C*temp2 + S*temp1);
            if((aura[3*i]-bodx)*(aura[3*i]-bodx) + (aura[3*i+1]-body)*(aura[3*i+1]-body)>0.1*0.1*0.7){
                auravel[2*i] = -aura[3*i] + bodx + 0.3*((static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5);
                auravel[2*i+1] = -aura[3*i+1]+body + 0.3*((static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5);
            }
        } else if(auradata[i] == UNBOUND){
            if ((mousex - aura[3*i])*(mousex - aura[3*i])
                + (mousey - aura[3*i+1])*(mousey - aura[3*i+1]) < 0.007){
                auradata[i] = DEAD;
            } else {
                temp1 = (mousex - aura[3*i])*(mousex - aura[3*i])
                + (mousey - aura[3*i+1])*(mousey - aura[3*i+1]);
                auravel[2*i] += 3*dt*(mousex - aura[3*i])/temp1;
                auravel[2*i+1] += 3*dt*(mousey - aura[3*i+1])/temp1;
                auravel[2*i] *= 1-0.005;
                auravel[2*i+1] *= 1-0.005;
            }
        } else {
            if(count > 0){
                auradata[i] = BOUND;
                aura[3*i] = bodx;
                aura[3*i+1] = body;
                auravel[2*i] = 0.3*((static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5);
                auravel[2*i+1] = 0.3*((static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5);
                count--;
            }
        }
    }
}

void MoveAura(float aura[], float dx, float dy, int auradata[]){
    for(int i = 0; i < auranum; i++){
        if(auradata[i] == BOUND){
            aura[3*i]+=dx;
            aura[3*i+1]+=dy;
        }
    }
}
void MoveUnbound(float aura[], float dx, float dy, int auradata[]){
    for(int i = 0; i < auranum; i++){
        if(auradata[i] == UNBOUND){
            aura[3*i]+=dx;
            aura[3*i+1]+=dy;
        }
    }
}

void Unbind(float aura[], int auradata[], float mousex, float mousey,
            float aurarot, float bodx, float body, float auravel[], float bodang, float bodspe){
    float min = 0;
    int mini = 0;
    bool first = true;
    float d;
    for(int i= 0; i < auranum; i++){
        if(auradata[i]==BOUND){
            d = (aura[3*i]-mousex)*(aura[3*i]-mousex) + (aura[3*i+1]-mousey)*(aura[3*i+1]-mousey);
            if( first || d < min ){
                first = false;
                min = d;
                mini = i;
            }
        }
    }
    auradata[mini] = UNBOUND;
    auravel[2*mini] -= aurarot*(aura[3*mini+1]-body);
    auravel[2*mini+1] += aurarot*(aura[3*mini]-bodx);
    auravel[2*mini] += bodspe* cos(bodang);
    auravel[2*mini+1] += bodspe* sin(bodang);
    
    
}

void MoveWalls(float Walls[], float dx, float dy, float Centers[], int NW, int NC){
    for(int i = 0; i < 3*NW; i++){
        Walls[2*i] += dx;
        Walls[2*i+1] += dy;
    }
    for(int i = 0; i < NC; i++){
        Centers[3*i] += dx;
        Centers[3*i+1] += dy;
    }
}

void MoveBullets(float bullets[], int bulletnum, float dx, float dy){
    for(int i = 0; i < bulletnum; i++){
        bullets[3*i] += dx;
        bullets[3*i + 1] += dy;
    }
}

void MoveEverything(float particles[], float body[], int NP, int NB, float dx, float dy,float background[], float density, float AspectRatio, float aura[], int auradata[], float Walls[],float Centers[], int NW, int NC, Monster1 m1[], int m1num,
                    float bullets[], int bulletnum){
    MoveParticles(particles, NP, dx, dy);
    MoveBody(body, NB, dx, dy);
    MoveBackground(background, dx, dy, density, AspectRatio);
    MoveAura(aura, dx, dy, auradata);
    MoveWalls(Walls, dx, dy, Centers, NW, NC);
    MoveBullets(bullets, bulletnum, dx,dy);
    for(int i = 0; i < m1num; i++){
        m1[i].move(dx,dy);
    }
}
void MoveRelative(float particles[], int NP, float dx, float dy,float background[],
                  float density, float AspectRatio, int auradata[],float aura[], float Walls[],float Centers[],
                  int NW, int NC, Monster1 m1[], int m1num, float bullets[], int bulletnum){
    MoveParticles(particles, NP, dx, dy);
    MoveBackground(background, dx, dy, density, AspectRatio);
    MoveUnbound( aura, dx,dy,auradata);
    MoveWalls(Walls, dx, dy, Centers, NW, NC);
    MoveBullets(bullets, bulletnum,dx,dy);
    for(int i = 0; i < m1num; i++){
        m1[i].move(dx,dy);
    }
}


int CollisionDetector(float aura[], float auravel[], int auradata[],
                      GLuint UnboundIndex[], int unb, float colnums[]
                      ,float bodx, float body, float aurarot, float vx, float vy){
    
    // might make a faster algorithm if lag becomes a problem
    
    int index = 0;
    
    bool UnDetected[auranum];
    
    for(int i = 0; i < auranum; i++){
        UnDetected[i] = true;
    }
    
    for(int i = 0; i < unb; i++){
        for(int j = 0; j < auranum; j++){
            float d = (aura[3*UnboundIndex[i]] - aura[3*j])*(aura[3*UnboundIndex[i]] - aura[3*j])
            + (aura[3*UnboundIndex[i]+1] - aura[3*j+1])*(aura[3*UnboundIndex[i]+1] - aura[3*j+1]);
            if(d < 0.0009765625 && UnDetected[j] && UnDetected[UnboundIndex[i]]){
                float dot =auravel[2*UnboundIndex[i]]*auravel[2*j]
                +auravel[2*UnboundIndex[i]+1]*auravel[2*j+1];
                if(auradata[j] == BOUND){
                    dot += aurarot*((aura[3*j]-bodx )*auravel[2*UnboundIndex[i]+1]
                                   +(body-(aura[3*j+1]) )*auravel[2*UnboundIndex[i]]);
                    dot += vy*auravel[2*UnboundIndex[i]+1] + vx*auravel[2*UnboundIndex[i]];
                }else if (auradata[j] == UNBOUND){
                    dot *= 10;
                }else {
                    dot = 1;
                }
                if(dot < -1){
                    UnDetected[j] = false;
                    UnDetected[UnboundIndex[i]] = false;
                    colnums[4* index] = 0.5 * (aura[3*UnboundIndex[i]] + aura[3*j]);
                    colnums[4* index + 1] = 0.5 * (aura[3*UnboundIndex[i] + 1] + aura[3*j + 1]);
                
                    colnums[4* index + 2] = -dot;
            
                    
                    if(auradata[j] == BOUND){
                        colnums[4* index + 3] = 0;
                    }else if (auradata[j] == UNBOUND){
                        colnums[4* index + 3] = 1;
                        auradata[j] = DEAD;
                    }else{
                        colnums[4* index + 3] = 2;
                    }
                    
                    auradata[UnboundIndex[i]] = DEAD;
                    index++;
                }
            }
        }
    }
    
    return index;
}

void AddBullet(float bullets[],float bulletvel[], int* bulletnum, int bulletmax, float x, float y, float vx, float vy){
    for(int i = *bulletnum; i > 0; i--){
        if(i < bulletmax){
            bullets[3*i] = bullets[3*(i-1)];
            bullets[3*i+1] = bullets[3*(i-1)+1];
            bullets[3*i+2] = bullets[3*(i-1)+2];
            bulletvel[2*i] = bulletvel[2*(i-1)];
            bulletvel[2*i+1] = bulletvel[2*(i-1)+1];
        }
    }
    bullets[0] = x;
    bullets[1] = y;
    bullets[2] = 0.01;
    bulletvel[0] = vx;
    bulletvel[1] = vy;
    if(*bulletnum < bulletmax){
        *bulletnum += 1;
    }
}

void RemoveBullet(float bullets[],float bulletvel[], int* bulletnum, int bulletmax, int index){
    for(int i = index; i < bulletmax-1; i++){
        bullets[3*i] = bullets[3*(i+1)];
        bullets[3*i+1] = bullets[3*(i+1)+1];
        bullets[3*i+2] = bullets[3*(i+1)+2];
        bulletvel[2*i] = bulletvel[2*(i+1)];
        bulletvel[2*i+1] = bulletvel[2*(i+1)+1];
    }
    if(*bulletnum > 0){
        *bulletnum -= 1;
    }
    
}

void IncrementBullets(float bullets[],float bulletvel[], int bulletnum, float dt){
    for(int i = 0; i < bulletnum; i++){
        bullets[3*i] += dt* bulletvel[2*i];
        bullets[3*i+1] += dt* bulletvel[2*i+1];
    }
}

void InflictMinor(float* HP, float* minorHP, float damage){
    if(damage< *HP){
        *HP -= damage;
        *minorHP += damage;
    }else if(damage < 2*(*HP) + *minorHP){
        *minorHP += *HP;
        *HP = 0;
        *minorHP -= damage;
    }else{
        *HP = 0;
        *minorHP = 0;
    }
}

void InflictMajor(float* HP, float* minorHP, float damage){
    if(damage < *minorHP){
        *minorHP -= damage;
    }else if (damage < 2*(*HP) + *minorHP){
        *minorHP = 0;
        *HP -= (damage - *minorHP)/2;
    }else{
        *HP = 0;
        *minorHP = 0;
    }
}


int main(void)
{
    
    float bodyrad = 0.02;
    
    // don't mess with the following -------------------------
    GLFWwindow* window;
    
    if (!glfwInit())
        return -1;
    
    const GLFWvidmode * mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    int W = mode->width;
    int H = mode->height;
    
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
     
    window = glfwCreateWindow(W, H, "Hello World", glfwGetPrimaryMonitor(), NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    

    glfwMakeContextCurrent(window);
    
    glfwSwapInterval(1);
    
    float AspectRatio = ((float) W)/((float) H);
    
    //-----------------------------------------------------------
    
    glEnable(GL_BLEND);
    
    // define vars
    float HP = 1;
    float minorHP = 0;
    
    
    float colnums[auranum*4];
    float explosions[auranum * 3];
    int explosiondata[auranum];
    int numexp;
    float particles[3* bodyparticles];
    float parvel[2* bodyparticles];
    float aura[3*auranum];
    float auravel[2*auranum];
    int auradata[auranum];
    float aurarot = 10;
    float randrad;
    float randang;
    for(int i = 0; i < 3* bodyparticles; i++){
        if(i%3 == 2){
            particles[i] = 0.01;
        }else{
            particles[i] = 0;
        }
    }
    for(int i = 0; i < auranum; i++){
        aura[i*3 + 2] = 0.015625;
        randrad = 0.1*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
        randang = 2* 3.1415926 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
        aura[i*3] = randrad *cos(randang);
        aura[i*3 + 1] = randrad *sin(randang) ;
    }
    for(int i = 0; i < auranum; i++){
        randrad = 0.3*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
        randang = 2* 3.1415926 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
        randrad += 0.01;
        auravel[i*2] = randrad *cos(randang) ;
        auravel[i*2 + 1] = randrad *sin(randang) ;
    }
    for(int i = 0; i < auranum; i++){
        auradata[i] = BOUND;
    }
    
    // define and bind vertex buffers
    
    GLuint vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    
    // Activate initial vertex buffer
    
    // glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * bodyparticles, particles, GL_DYNAMIC_DRAW);
    
    
    // define vertex shaders
    
    static const char* vertex_shader_text =
    "#version 410 core\n"
    "uniform float AR;"
    "in vec2 position;"
    "in float radius;"
    "out float Radius;"
    "void main()"
    "{"
    "    Radius = radius;"
    "    gl_Position = vec4(position.x, position.y * AR , 0.0, 1.0);"
    "}";
    
    static const char* vertex_shader_text_2 =
    "#version 410 core\n"
    "in vec2 position;"
    "in vec2 gradient;"
    "out vec2 Grad;"
    "void main()"
    "{"
    "    Grad = gradient;"
    "    gl_Position = vec4(position, 0.0, 1.0);"
    "}";
    
    static const char* vertex_shader_text_3 =
    "#version 410 core\n"
    "uniform float AR;"
    "in vec2 position;"
    "void main()"
    "{"
    "    gl_Position = vec4(position.x, AR* position.y, 0.0, 1.0);"
    "}";
    
    
    // define geometry shaders
     
    static const char* geometry_shader_text =
    "#version 410 core\n"
    "const float PI = 3.1415926;"
    "uniform float AR;"
    "layout(points) in;"
    "layout(triangle_strip, max_vertices = 20) out;"
    "in float Radius[];"
    "void main()"
    "{"
    "for(int i = 0; i < 10; i++){"
    "gl_Position = gl_in[0].gl_Position + Radius[0] * vec4(sin((PI / 10) * i), AR * cos((PI / 10) * i), 0.0, 0.0);"
    "EmitVertex();"
    "gl_Position = gl_in[0].gl_Position + Radius[0] * vec4(sin ( - (PI / 10) * ( i + 1 ) ), AR * cos((PI / 10) * ( i + 1 )), 0.0, 0.0);"
    "EmitVertex();"
    "}"
    "EndPrimitive();"
    "}";
    
    static const char* geometry_shader_text_2 =
    "#version 410 core\n"
    "float smthstp(float x){"
    "return x*x*(3 - 2*x);"
    "}"
    "layout(triangles) in;"
    "layout(triangle_strip, max_vertices = 24) out;"
    "in vec2 Grad[];"
    "out float intensity;"
    "void main()"
    "{"
    "float X;"
    "float Y;"
    "float Z;"
    "vec4 A = gl_in[0].gl_Position;"
    "vec4 B = gl_in[1].gl_Position;"
    "vec4 C = gl_in[2].gl_Position;"
    "for(int i = 1; i < 5; i = i+1){"
    "for(int j = 0; j < 2*i + 1; j = j+1){"
    "float a = float(i);"
    "float b = float(j);"
    "X = ((2*a - b) - mod( b , 2))/8;"
    "Y = (b - mod(b , 2)) /8;"
    "Z = 1 - ( a  -  mod( b , 2) )/4;"
    "vec4 pos = X * A  + Y * B  + Z * C;"
    "gl_Position = pos;"
    "intensity = smthstp(X) * dot((pos - A).xy,Grad[0]) + smthstp(Y) * dot((pos - B).xy,Grad[1]) + smthstp(Z) * dot((pos - C).xy,Grad[2]);"
    "EmitVertex();"
    "}"
    "EndPrimitive();"
    "}"
    "}";
    
    static const char* geometry_shader_text_3 =
    "#version 410 core\n"
    "const float PI = 3.1415926;"
    "uniform float sides;"
    "layout(triangles) in;"
    "layout(triangle_strip, max_vertices = 300) out;"
    "in float Radius[];"
    "void main()"
    "{"
    "for(int i = 0; i < sides; i++){"
    "gl_Position = gl_in[0].gl_Position + gl_in[2].gl_Position - gl_in[1].gl_Position + (gl_in[1].gl_Position - gl_in[0].gl_Position) * sin((PI / (2* sides)) * i) + (gl_in[1].gl_Position - gl_in[2].gl_Position) * cos((PI / (2* sides)) * i);"
    "EmitVertex();"
    "gl_Position = gl_in[0].gl_Position + gl_in[2].gl_Position - gl_in[1].gl_Position + (gl_in[1].gl_Position - gl_in[0].gl_Position) * sin((PI / (2* sides)) * (i+1)) + (gl_in[1].gl_Position - gl_in[2].gl_Position) * cos((PI / (2* sides)) * (i+1));"
    "EmitVertex();"
    "gl_Position = gl_in[1].gl_Position;"
    "EmitVertex();"
    "EndPrimitive();"
    "}"
    "}";
    
    
    // define fragment shaders
    
    static const char* fragment_shader_text =
    "#version 410 core\n"
    "uniform vec3 triangleColor;"
    "uniform float opacity;"
    "out vec4 outColor;"
    "void main()"
    "{"
    "    outColor = vec4(triangleColor, opacity);"
    "}";
    
    static const char* fragment_shader_text_2 =
    "#version 410 core\n"
    "uniform vec3 Color0;"
    "uniform vec3 Color1;"
    "uniform float scale;"
    "in float intensity;"
    "out vec4 outColor;"
    "void main()"
    "{"
    "    outColor = vec4( Color0 + ((atan(scale * intensity) / 3.1415926) + 0.5) * Color1  , 1.0);"
    "}";
    
    static const char* fragment_shader_text_3 =
    "#version 410 core\n"
    "uniform vec3 triangleColor;"
    "uniform float opacity;"
    "out vec4 outColor;"
    "void main()"
    "{"
    "    outColor = vec4(triangleColor, opacity);"
    "}";
    
    
    // compile the shaders and check errors
    
    GLint status;
    
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertex_shader_text, NULL);
    glCompileShader(vertexShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled vertex shader correctly" << endl;
    }else{
        cout << "Vertex shader compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(vertexShader, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    GLuint vertexShader2 = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader2, 1, &vertex_shader_text_2, NULL);
    glCompileShader(vertexShader2);
    glGetShaderiv(vertexShader2, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled vertex shader 2 correctly" << endl;
    }else{
        cout << "Vertex shader 2 compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(vertexShader2, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    GLuint vertexShader3 = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader3, 1, &vertex_shader_text_3, NULL);
    glCompileShader(vertexShader3);
    glGetShaderiv(vertexShader3, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled vertex shader 3 correctly" << endl;
    }else{
        cout << "Vertex shader 3 compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(vertexShader3, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    
    GLuint geometryShader = glCreateShader(GL_GEOMETRY_SHADER);
    glShaderSource(geometryShader, 1, &geometry_shader_text, NULL);
    glCompileShader(geometryShader);
    glGetShaderiv(geometryShader, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled geometry shader correctly" << endl;
    }else{
        cout << "Geometry shader compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(geometryShader, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    GLuint geometryShader2 = glCreateShader(GL_GEOMETRY_SHADER);
    glShaderSource(geometryShader2, 1, &geometry_shader_text_2, NULL);
    glCompileShader(geometryShader2);
    glGetShaderiv(geometryShader2, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled geometry shader 2 correctly" << endl;
    }else{
        cout << "Geometry shader 2 compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(geometryShader2, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    GLuint geometryShader3 = glCreateShader(GL_GEOMETRY_SHADER);
    glShaderSource(geometryShader3, 1, &geometry_shader_text_3, NULL);
    glCompileShader(geometryShader3);
    glGetShaderiv(geometryShader3, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled geometry shader 3 correctly" << endl;
    }else{
        cout << "Geometry shader 3 compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(geometryShader3, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragment_shader_text, NULL);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled fragment shader correctly" << endl;
    }else{
        cout << "Fragment shader compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(fragmentShader, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    GLuint fragmentShader2 = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader2, 1, &fragment_shader_text_2, NULL);
    glCompileShader(fragmentShader2);
    glGetShaderiv(fragmentShader2, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled fragment shader 2 correctly" << endl;
    }else{
        cout << "Fragment shader 2 compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(fragmentShader2, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    GLuint fragmentShader3 = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader3, 1, &fragment_shader_text_3, NULL);
    glCompileShader(fragmentShader3);
    glGetShaderiv(fragmentShader3, GL_COMPILE_STATUS, &status);
    if (status == GL_TRUE){
        cout << "Compiled fragment shader 3 correctly" << endl;
    }else{
        cout << "Fragment shader 3 compilation error" << endl;
        char buffer[512];
        glGetShaderInfoLog(fragmentShader3, 512, NULL, buffer);
        
        for(int i = 0; i < 512; i++){
            cout << buffer[i];
        }
        cout<< endl;
        cout<< "OpenGL version:   ";
        cout << glGetString( GL_VERSION ) << endl;
    }
    
    
    
    // Create programs
    
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, geometryShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram); // must be re-called if shader is modified.
    
    GLuint shaderProgram2 = glCreateProgram();
    glAttachShader(shaderProgram2, vertexShader2);
    glAttachShader(shaderProgram2, geometryShader2);
    glAttachShader(shaderProgram2, fragmentShader2);
    glLinkProgram(shaderProgram2);
    
    GLuint shaderProgram3 = glCreateProgram();
    glAttachShader(shaderProgram3, vertexShader3);
    glAttachShader(shaderProgram3, geometryShader3);
    glAttachShader(shaderProgram3, fragmentShader3);
    glLinkProgram(shaderProgram3);
    
    
    // Enable the attributes
    
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    GLuint ebo;
    glGenBuffers(1, &ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, 3*sizeof(float), 0);
    
    GLint radAttrib = glGetAttribLocation(shaderProgram, "radius");
    glEnableVertexAttribArray(radAttrib);
    glVertexAttribPointer(radAttrib, 1, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)(2*sizeof(float)));
    
    GLuint vao2;
    glGenVertexArrays(1, &vao2);
    glBindVertexArray(vao2);
    
    GLuint ebo2;
    glGenBuffers(1, &ebo2);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo2);
    
    GLint posAttrib2 = glGetAttribLocation(shaderProgram2, "position");
    glEnableVertexAttribArray(posAttrib2);
    glVertexAttribPointer(posAttrib2, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), 0);
    
    GLint gradAttrib = glGetAttribLocation(shaderProgram2, "gradient");
    glEnableVertexAttribArray(gradAttrib);
    glVertexAttribPointer(gradAttrib, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)(2*sizeof(float)));
    
    GLuint vao3;
    glGenVertexArrays(1, &vao3);
    glBindVertexArray(vao3);
    
    GLint posAttrib3 = glGetAttribLocation(shaderProgram3, "position");
    glEnableVertexAttribArray(posAttrib3);
    glVertexAttribPointer(posAttrib3, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), 0);
    
    
    // Uniform variables unf
    
    glUseProgram(shaderProgram);
    GLint uniColor = glGetUniformLocation(shaderProgram, "triangleColor");
    GLint opacity = glGetUniformLocation(shaderProgram, "opacity");
    GLint ARuniform = glGetUniformLocation(shaderProgram, "AR");
    glUniform3f(uniColor, 1.0f, 1.0f, 1.0f);
    glUniform1f(ARuniform, AspectRatio);
    glUniform1f(opacity, 1);
    glUseProgram(shaderProgram2);
    GLint col0 = glGetUniformLocation(shaderProgram2, "Color0");
    GLint col1 = glGetUniformLocation(shaderProgram2, "Color1");
    GLint scale = glGetUniformLocation(shaderProgram2, "scale");
    glUniform3f(col0, 0.185,0,0);
    glUniform3f(col1, 0.7,1.0,1.0);
    glUniform1f(scale, 3);
    glUseProgram(shaderProgram3);
    GLint uniColor3 = glGetUniformLocation(shaderProgram3, "triangleColor");
    GLint opacity3 = glGetUniformLocation(shaderProgram3, "opacity");
    GLint ARuniform3 = glGetUniformLocation(shaderProgram3, "AR");
    GLint sides3 = glGetUniformLocation(shaderProgram3, "sides");
    glUniform3f(uniColor3, 0.1, 0.1, 0.1);
    glUniform1f(ARuniform3, AspectRatio);
    glUniform1f(opacity3, 1);
    glUniform1f(sides3, 20);
    

    // vars
    
    GLuint UnboundIndex[auranum];
    GLuint BoundIndex[auranum];
    
    float regentime = 0;
    float prevtime20 = 0;
    float prevtime = 0;
    
    float body[60];
    for(int i = 0; i < 20; i++){
        body[3*i] = 0;
        body[3*i + 1] = 0;
        body[3*i + 2] = bodyrad * i / 20;
    }
    
    float bodang = 0;
    float bodspe = 0;
    
    int w;
    int a;
    int s;
    int d;
    int z;
    int x;
    int click;
    int space;
    w = glfwGetKey(window, GLFW_KEY_W);
    a = glfwGetKey(window, GLFW_KEY_A);
    s = glfwGetKey(window, GLFW_KEY_S);
    d = glfwGetKey(window, GLFW_KEY_D);
    z = glfwGetKey(window, GLFW_KEY_Z);
    x = glfwGetKey(window, GLFW_KEY_X);
    click = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    space = glfwGetKey(window, GLFW_KEY_SPACE);
    
    float maxspeed = 1;
    
    float Wacc = 0.5;
    
    float ADvel = 4;
    
    float Sdecel = 0.9;
    
    float framedx = 0;
    float framedy = 0;
    
    double mousex = 0;
    double mousey = 0;
    
    float background[4*9*9];
    
    float density = 1.7320508;
    
    for(int i = 0; i < 9; i++){
        for(int j = 0; j < 9; j++){
            background[4*9*i + 4*j] = (((float)(i-4))/density) * 0.8660254;
            background[4*9*i + 4*j + 1] = (((float)(j-4)) + 0.5 * ((float) (i % 2)))* AspectRatio /density;
            
            background[4*9*i + 4*j + 2] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
            background[4*9*i + 4*j + 3] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
        }
    }
    

    GLuint elements[2*8*8*3];
    
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            elements[16*3*i + 3*j] = 9*i + j;
            elements[16*3*i + 3*j + 1] = 9*i + j + 1;
            elements[16*3*i + 3*j + 2] = 9*i + j + 9;
        }
    }
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            elements[8*3*(2*i+1) + 3*j] = 9*(i+1) + j;
            elements[8*3*(2*i+1) + 3*j + 1] = 9*(i+1) + j + 1;
            elements[8*3*(2*i+1) + 3*j + 2] = 9*(i+1) + j - 8;
        }
    }
    
    float C = cos(0.05);
    float S = sin(0.05);
    numexp = 0;
    
    int regen = 0;
    
    // LEVEL DATA:
    
    float Walls[] = {-1,0,-1,-1,0,-1, 0,-1,1,-1,1,0, 0,1,-1,1,-1,0, 1,0,2,0,2,1, 2,1,2,2,1,2, 1,2,0,2,0,1,
        -1,1,0,1,0,2, 1,-1,1,0,2,0, -1,1,0,2,-1,2, 1,-1,2,0,2,-1,
        -1,-1,-8,-1,-1,6, 2,-1,2,-8,-5,-1, 2,2,9,2,2,-5, -1,2,-1,9,6,2};
    float Centers[] = {0,0,1,1,1,1};
    float WallsCopy[] = {-1,0,-1,-1,0,-1, 0,-1,1,-1,1,0, 0,1,-1,1,-1,0, 1,0,2,0,2,1, 2,1,2,2,1,2, 1,2,0,2,0,1,
        -1,1,0,1,0,2, 1,-1,1,0,2,0, -1,1,0,2,-1,2, 1,-1,2,0,2,-1,
        -1,-1,-8,-1,-1,6, 2,-1,2,-8,-5,-1, 2,2,9,2,2,-5, -1,2,-1,9,6,2};
    float CentersCopy[] = {0,0,1,1,1,1};
    int wallnum = 14;
    int centnum = 2;
    int CurveStraightPartition = 6;
    
    float Fade = 0;
    
    Monster1 m1[2];
    float bullets[300];
    float bulletvel[200];
    for(int i = 0; i < 100; i++){
        bullets[3*i] = 0;
        bullets[3*i+1] = 0;
        bullets[3*i+2] = 0;
        bulletvel[2*i] = 0;
        bulletvel[2*i+1] = 0;
    }
    int bulletnum=0;
    int m1num =  2;
    m1[0].centerx0 = 1-0.2;
    m1[0].centery0 = 1;
    m1[0].territoryRad = 1-0.2;
    m1[1].centerx0 = 1+0.2;
    m1[1].centery0 = 1;
    m1[1].territoryRad = 1-0.2;
    m1[0].respawn();
    m1[1].respawn();
    
    
    // MAIN LOOP ___________________
    
    
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);
        
        
        // --------- drawing and incrementation goes here --------
        
        
        float time = (float) glfwGetTime();
        
        // regen loop
        
        regen = 0;
        while(time - regentime > 0.1){
            if(minorHP >= 0.002){
                HP+= 0.002;
                minorHP-= 0.002;
            }else if(minorHP>0){
                HP += minorHP;
                minorHP = 0;
            }
            for(int mn = 0; mn < m1num; mn++){
                if(m1[mn].minorHP >= 0.002){
                    m1[mn].HP+= 0.002;
                    m1[mn].minorHP-= 0.002;
                }else if(m1[mn].minorHP>0){
                    m1[mn].HP += m1[mn].minorHP;
                    m1[mn].minorHP = 0;
                }
            }
            regen += 1;
            regentime = time;
        }
        
        
        // 20 times per second loop:
        
        while(time - prevtime20 > 0.05){
            prevtime20 += 0.05;
            IncrementParticles(particles,parvel,bodyparticles,0.05,body[0],body[1],
                               bodang + 3.1415926,0.012,0.065,0.2);
            IncrementBackground(background,C,S);
            if (Fade > 0) {
                Fade -= 0.01;
            }
            int numexpcopy = numexp;
            for(int i = 0; i < numexpcopy; i++){
                explosions[3*i+2] += 0.03;
                if(explosions[3*i+2] > 0.15){
                    numexp--;
                }
            }
            
            w = glfwGetKey(window, GLFW_KEY_W);
            a = glfwGetKey(window, GLFW_KEY_A);
            s = glfwGetKey(window, GLFW_KEY_S);
            d = glfwGetKey(window, GLFW_KEY_D);
            z = glfwGetKey(window, GLFW_KEY_Z);
            x = glfwGetKey(window, GLFW_KEY_X);
            space = glfwGetKey(window, GLFW_KEY_SPACE);
            if (w == GLFW_PRESS && bodspe < maxspeed){
                bodspe += Wacc * 0.05;
            }
            if (a == GLFW_PRESS){
                bodang += ADvel * 0.05;
            }
            if (d == GLFW_PRESS){
                bodang -= ADvel * 0.05;
            }
            if (s == GLFW_PRESS){
                bodspe *= Sdecel;
            }
            if (z == GLFW_PRESS && aurarot <19){
                aurarot += 3;
            }
            if (x == GLFW_PRESS && aurarot > -19){
                aurarot -= 3;
            }
            
            click = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
            if (click == GLFW_PRESS){
                Unbind(aura, auradata,(2*((float)mousex)/W - 1),-(2*((float)mousey)/H - 1),
                       aurarot,body[0],body[1],auravel,bodang,bodspe);
            }
        }
        
        // Bounce against walls
        
        bool OutOfBounds = true;
        int NN = 0;
        float min = 0;
        bool first = true;
        for(int i = 0; i < centnum; i++){
            if((Centers[3*i]-body[0])*(Centers[3*i]-body[0]) + (Centers[3*i+1]-body[1])*(Centers[3*i+1]-body[1])
               > Centers[3*i+2]*Centers[3*i+2]){
                if(((Centers[3*i]-body[0])*(Centers[3*i]-body[0]) + (Centers[3*i+1]-body[1])*(Centers[3*i+1]-body[1]) < min || first)
                   && sin(bodang)*(Centers[3*i+1]-body[1]) + cos(bodang)*(Centers[3*i]-body[0]) < 0){
                    min = (Centers[3*i]-body[0])*(Centers[3*i]-body[0]) + (Centers[3*i+1]-body[1])*(Centers[3*i+1]-body[1]);
                    NN = i;
                    first = false;
                }
            }else{
                OutOfBounds = false;
            }
        }
        
        if(OutOfBounds){
            if(bodspe >= 0.9){
                InflictMajor(&HP, &minorHP, 0.01);
            }
            InflictMinor(&HP, &minorHP, 0.03);
            bodang = 2*atan2(-(Centers[3*NN]-body[0]),(Centers[3*NN+1]-body[1])) - bodang;
            MoveRelative(particles, bodyparticles,
                         -(Centers[3*NN]-body[0])*0.01, -(Centers[3*NN+1]-body[1])*0.01,
                         background, density, AspectRatio, auradata, aura, Walls, Centers,
                         wallnum, centnum, m1, m1num, bullets, bulletnum);
        }
        
        
        
        
        // frame by frame motion
        
        glfwGetCursorPos(window, &mousex, &mousey);
        framedx = ((2*((float)mousex)/W - 1) + body[0])/2;
        framedy = -((2*((float)mousey)/H - 1)/AspectRatio - body[1])/2;
        IncrementAura(aura, auravel, auradata, aurarot,(time - prevtime),body[0], body[1],
                    (2*((float)mousex)/W - 1),-(2*((float)mousey)/H - 1)/AspectRatio, regen);
        MoveEverything(particles,body,bodyparticles,20,-framedx,-framedy,
                       background,density,AspectRatio,aura,auradata, Walls, Centers, wallnum,
                       centnum,m1,m1num,bullets,bulletnum);
        MoveRelative(particles, bodyparticles,
                     -bodspe * cos(bodang) * (time - prevtime), - bodspe * sin(bodang) * (time - prevtime),
                     background, density, AspectRatio, auradata, aura, Walls, Centers, wallnum,
                     centnum,m1,m1num,bullets,bulletnum);
        
        for(int i = 0; i < m1num; i++){
            if(m1[i].alive){
                m1[i].increment(body[0], body[1], bodspe * cos(bodang), bodspe * sin(bodang), time - prevtime);
                for(int j = 0; j < 5; j++){
                    if(m1[i].recharge[j] == 0){
                        AddBullet(bullets, bulletvel, &bulletnum, 100, m1[i].x + 0.07*cos(m1[i].sang + (2*3.1415926/5)*j), m1[i].y + 0.07*sin(m1[i].sang + (2*3.1415926/5)*j), 1.5*cos(m1[i].sang + (2*3.1415926/5)*j), 1.5*sin(m1[i].sang + (2*3.1415926/5)*j));
                    }
                }
                for(int j = 0; j < numexp; j++){
                    if((explosions[3*j]-m1[i].x)*(explosions[3*j]-m1[i].x) + (explosions[3*j+1]-m1[i].y)*(explosions[3*j+1]-m1[i].y)
                       < (explosions[3*j+2]+0.07)*(explosions[3*j+2]+0.07)){
                        if(explosiondata[j]==0){
                            InflictMajor(&m1[i].HP, &m1[i].minorHP, (time - prevtime));
                            InflictMinor(&m1[i].HP, &m1[i].minorHP, (time - prevtime)*2);
                        }
                        if(explosiondata[j] == 1){
                            InflictMinor(&m1[i].HP, &m1[i].minorHP, (time - prevtime)/2);
                        }
                    }
                }
            }
        }
        
        for(int j = 0; j < bulletnum; j++){
            OutOfBounds = true;
            for(int i = 0; i < centnum; i++){
                if((Centers[3*i]-bullets[3*j])*(Centers[3*i]-bullets[3*j]) + (Centers[3*i+1]-bullets[3*j+1])*(Centers[3*i+1]-bullets[3*j+1])
                   < Centers[3*i+2]*Centers[3*i+2]){
                    OutOfBounds = false;
                }
            }
            if(OutOfBounds){
                RemoveBullet(bullets, bulletvel, &bulletnum, 100, j);
            }
        }
        for(int j = 0; j < bulletnum; j++){
            if((body[0]-bullets[3*j])*(body[0]-bullets[3*j]) + (body[1]-bullets[3*j+1])*(body[1]-bullets[3*j+1])
               < (bodyrad + 0.01)*(bodyrad + 0.01)){
                RemoveBullet(bullets, bulletvel, &bulletnum, 100, j);
                InflictMajor(&HP, &minorHP, 0.05);
                InflictMinor(&HP, &minorHP, 0.05);
                
            }
        }
        
        
        
        
        
        IncrementBullets(bullets, bulletvel, bulletnum, time-prevtime);
        prevtime = time;
        
        
        
        
        // Collision Detection follows
        int bnd = 0;
        int unb = 0;
        for(int i = 0; i < auranum; i++){
            if(auradata[i] == BOUND){
                BoundIndex[bnd] = i;
                bnd++;
            }
            if(auradata[i] == UNBOUND){
                UnboundIndex[unb] = i;
                unb++;
            }
        }
        
        
        if (space == GLFW_PRESS){
            int newexp;
            newexp = CollisionDetector(aura, auravel, auradata, UnboundIndex, unb, colnums,
                                       body[0], body[1], 2 * aurarot , bodspe * cos(bodang),  bodspe * sin(bodang));
            for(int i = numexp-1; i >= 0; i--){
                if(i + newexp < auranum){
                    explosions[(i + newexp)*3] = explosions[i*3];
                    explosions[(i + newexp)*3+1] = explosions[i*3+1];
                    explosions[(i + newexp)*3+2] = explosions[i*3+2];
                    explosiondata[i + newexp] = explosiondata[i];
                }
            }
            for(int i = 0; i < newexp; i++){
                explosions[i*3] = colnums[i*4];
                explosions[i*3+1] = colnums[i*4+1];
                explosions[i*3+2] = 0.001;
                explosiondata[i] = (int)colnums[i*4+3];
            }
            numexp+= newexp;
        }
        
        // Draw everything now
        
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // background
        
        glUseProgram(shaderProgram2);
        glBindVertexArray(vao2);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4*9*9, background, GL_STATIC_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 2*8*8*3, elements, GL_STATIC_DRAW);
        glDrawElements(GL_TRIANGLES, 2*8*8*3, GL_UNSIGNED_INT, nullptr);
        
        
        
        glUseProgram(shaderProgram3);
        glBindVertexArray(vao3);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * wallnum , Walls, GL_STATIC_DRAW);
        glUniform3f(uniColor3, 0.1, 0.1, 0.1);
        glUniform1f(opacity3, 1);
        glUniform1f(sides3, 20);
        glDrawArrays(GL_TRIANGLES, 0, 3*CurveStraightPartition);
        glUniform3f(uniColor3, 0.1, 0.1, 0.1);
        glUniform1f(sides3, 1);
        glDrawArrays(GL_TRIANGLES,  3*CurveStraightPartition, 3*wallnum);
        
        // draw cuties
        for(int mn = 0; mn < m1num; mn++){
            if(m1[mn].alive){
                
                glUseProgram(shaderProgram);
                glBindVertexArray(vao);
                float Cutie[3] = {m1[mn].x,m1[mn].y,0.03};
                glUniform3f(uniColor, 0, 0.2, 0);  // color of circle in middle behind main body orig. 0.2,0.2,0.2
                glUniform1f(opacity, 1);
                glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 , Cutie, GL_DYNAMIC_DRAW);
                glDrawArrays(GL_POINTS, 0, 1 );
                
                glUseProgram(shaderProgram3);
                glBindVertexArray(vao3);
                float wing[6];
                glUniform1f(opacity3, 0.8);
                glUniform1f(sides3, 20);
                glUniform3f(uniColor3, 0.4, 0.3, 0.2); // color of main body orig. 0.65,0.22,0.29
                for(int i = 0; i < 5; i++){
                    wing[0] = m1[mn].x + 0.07* cos(m1[mn].sang + (2*3.1415926 / 5)*i);
                    wing[1] = m1[mn].y + 0.07* sin(m1[mn].sang + (2*3.1415926 / 5)*i);
                    wing[2] = m1[mn].x;
                    wing[3] = m1[mn].y;
                    wing[4] = m1[mn].x + 0.07* cos(m1[mn].sang + (2*3.1415926 / 5)*(i+1));
                    wing[5] = m1[mn].y + 0.07* sin(m1[mn].sang + (2*3.1415926 / 5)*(i+1));
                    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 , wing, GL_DYNAMIC_DRAW);
                    glDrawArrays(GL_TRIANGLES, 0, 3);
                }
                
                glUseProgram(shaderProgram);
                glBindVertexArray(vao);
                glUniform3f(uniColor, 0.3, 0, 0.2); // bullet color, orig. 0.3,0,0.2
                glUniform1f(opacity, 1);
                for(int i = 0; i < 5; i++){
                    Cutie[0] = m1[mn].x + 0.07* cos(m1[mn].sang + (2*3.1415926 / 5)*i);
                    Cutie[1] = m1[mn].y + 0.07* sin(m1[mn].sang + (2*3.1415926 / 5)*i);
                    Cutie[2] = 0.01 * m1[mn].recharge[i];
                    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 , Cutie, GL_DYNAMIC_DRAW);
                    glDrawArrays(GL_POINTS, 0, 1 );
                }
                glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3* bulletnum , bullets, GL_DYNAMIC_DRAW);
                glDrawArrays(GL_POINTS, 0, bulletnum);
                
                
                float LLCx = m1[mn].x -0.07;
                float LLCy = m1[mn].y + 0.1;
                float hght = 0.01;
                float lnth = 0.14;
                
                float GHE[] = {LLCx,LLCy,LLCx,LLCy+hght,LLCx+lnth*m1[mn].HP,LLCy,
                    LLCx+lnth*m1[mn].HP,LLCy,LLCx+lnth*m1[mn].HP,LLCy+hght,LLCx,LLCy+hght,};
                float YHE[] = {LLCx+lnth*m1[mn].HP,LLCy,LLCx+lnth*m1[mn].HP,LLCy+hght,LLCx+lnth*(m1[mn].HP+m1[mn].minorHP),LLCy,
                    LLCx+lnth*(m1[mn].HP+m1[mn].minorHP),LLCy,LLCx+lnth*(m1[mn].HP+m1[mn].minorHP),LLCy+hght,LLCx+lnth*m1[mn].HP,LLCy+hght};
                float RHE[] = {LLCx+lnth*(m1[mn].HP+m1[mn].minorHP),LLCy,LLCx+lnth*(m1[mn].HP+m1[mn].minorHP),LLCy+hght,LLCx+lnth,LLCy,
                    LLCx+lnth,LLCy,LLCx+lnth,LLCy+hght,LLCx+lnth*(m1[mn].HP+m1[mn].minorHP),LLCy+hght,};
                
                glUseProgram(shaderProgram3);
                glBindVertexArray(vao3);
                glUniform1f(sides3, 1);
                glUniform1f(opacity3, 0.7);
                glUniform3f(uniColor3, 0.3f, 0.8f, 0.3f);
                glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12 , GHE, GL_DYNAMIC_DRAW);
                glDrawArrays(GL_TRIANGLES,  0, 6);
                glUniform3f(uniColor3, 0.8f, 0.8f, 0.3f);
                glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12  , YHE, GL_DYNAMIC_DRAW);
                glDrawArrays(GL_TRIANGLES,  0, 6);
                glUniform3f(uniColor3, 0.8f, 0.3f, 0.3f);
                glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12  , RHE, GL_DYNAMIC_DRAW);
                glDrawArrays(GL_TRIANGLES,  0, 6);
            }
        }
        
        
        // circles
        glUseProgram(shaderProgram);
        glBindVertexArray(vao);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        
        
        bnd = 0;
        unb = 0;
        for(int i = 0; i < auranum; i++){
            if(auradata[i] == BOUND){
                BoundIndex[bnd] = i;
                bnd++;
            }
            if(auradata[i] == UNBOUND){
                UnboundIndex[unb] = i;
                unb++;
            }
        }
        
        glUniform3f(uniColor, 0.9f, 1.0f, 0.9f);
        glUniform1f(opacity, 0.15);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * auranum, aura, GL_DYNAMIC_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * auranum, UnboundIndex, GL_DYNAMIC_DRAW);
        glDrawElements(GL_POINTS, unb, GL_UNSIGNED_INT, nullptr);
        glUniform3f(uniColor, 1.0f, 1.0f, 1.0f);
        glUniform1f(opacity, 0.05);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * auranum, BoundIndex, GL_DYNAMIC_DRAW);
        glDrawElements(GL_POINTS, bnd, GL_UNSIGNED_INT, nullptr);
        
        
        glUniform3f(uniColor, 0.8f, 0.2f, 0.2f);
        glUniform1f(opacity, 0.2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * numexp , explosions, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_POINTS, 0, numexp );
        
         
        glUniform3f(uniColor, 1.0f, 1.0f, 1.0f);
        glUniform1f(opacity, 0.3);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * bodyparticles, particles, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_POINTS, 0,  bodyparticles);
        
        glUniform3f(uniColor, 1.0f, 1.0f, 1.0f);
        glUniform1f(opacity, 0.2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * 20 , body, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_POINTS, 0, 20 );
        
        float GH[] = {-1,1/AspectRatio - 0.005f, -1,1/AspectRatio - 0.02f,HP*0.99f-1,1/AspectRatio - 0.005f,
                    HP*0.99f-1,1/AspectRatio -0.005f,HP*0.99f-1,1/AspectRatio - 0.02f,-1,1/AspectRatio - 0.02f  };
        float YH[] = {HP*0.99f-1,1/AspectRatio - 0.005f, HP*0.99f-1,1/AspectRatio - 0.02f,HP*0.99f + minorHP*0.99f-1,1/AspectRatio - 0.005f,
            HP*0.99f + minorHP*0.99f-1,1/AspectRatio -0.005f,HP*0.99f + minorHP*0.99f-1,1/AspectRatio - 0.02f,HP*0.99f-1,1/AspectRatio - 0.02f  };
        float RH[] = {HP*0.99f + minorHP*0.99f-1,1/AspectRatio - 0.005f, HP*0.99f + minorHP*0.99f-1,1/AspectRatio - 0.02f,-0.01,1/AspectRatio - 0.005f,
            -0.01,1/AspectRatio -0.005f,-0.01,1/AspectRatio - 0.02f,HP*0.99f + minorHP*0.99f-1,1/AspectRatio - 0.02f  };
        
        glUseProgram(shaderProgram3);
        glBindVertexArray(vao3);
        glUniform1f(sides3, 1);
        glUniform1f(opacity3, 0.7);
        glUniform3f(uniColor3, 0.3f, 0.8f, 0.3f);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12 , GH, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES,  0, 6);
        glUniform3f(uniColor3, 0.8f, 0.8f, 0.3f);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12  , YH, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES,  0, 6);
        glUniform3f(uniColor3, 0.8f, 0.3f, 0.3f);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12  , RH, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES,  0, 6);
        
        float BB = ((float)bnd)/auranum;
        float UU = ((float)unb)/auranum;

        float BA[] = {1,1/AspectRatio - 0.005f, 1,1/AspectRatio - 0.02f,-BB*0.99f+1,1/AspectRatio - 0.005f,
            -BB*0.99f+1,1/AspectRatio -0.005f,-BB*0.99f+1,1/AspectRatio - 0.02f,1,1/AspectRatio - 0.02f  };
        float UA[] = {-BB*0.99f+1,1/AspectRatio - 0.005f, -BB*0.99f+1,1/AspectRatio - 0.02f,-BB*0.99f - UU*0.99f+1,1/AspectRatio - 0.005f,
            -BB*0.99f -UU*0.99f+1,1/AspectRatio -0.005f,-BB*0.99f - UU*0.99f+1,1/AspectRatio - 0.02f,-BB*0.99f+1,1/AspectRatio - 0.02f  };
        float DA[] = {-BB*0.99f -UU*0.99f+1,1/AspectRatio - 0.005f, -BB*0.99f + -UU*0.99f+1,1/AspectRatio - 0.02f,0.01,1/AspectRatio - 0.005f,
            0.01,1/AspectRatio -0.005f,0.01,1/AspectRatio - 0.02f,-BB*0.99f -UU*0.99f+1,1/AspectRatio - 0.02f  };
        
        
        glUniform1f(opacity3, 0.7);
        glUniform3f(uniColor3, 0.5f, 0.5f, 1.0f);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12, BA, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES,  0, 6);
        glUniform3f(uniColor3, 0.2f, 0.2f, 0.8f);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12 , UA, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES,  0, 6);
        glUniform3f(uniColor3, 0, 0, 0);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12 , DA, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES,  0, 6);
        
        
        float deathscreen[] = {-1,-1,-1,1,1,-1,-1,1,1,-1,1,1};
        glUniform3f(uniColor3, 0.3, 0, 0);
        glUniform1f(opacity3, Fade);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12 , deathscreen, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES,  0, 6);
        

        
        
        if(HP <= 0 &&  minorHP <= 0){
            for(int i = 0; i < 6*wallnum; i++){
                Walls[i] = WallsCopy[i];
            }
            for(int i = 0; i < 3*centnum; i++){
                Centers[i] = CentersCopy[i];
            }
            HP = 1;
            minorHP = 0;
            for(int i = 0; i < 20; i++){
                body[3*i] = 0;
                body[3*i+1] = 0;
            }
            bodspe = 0;
            bodang = 0;
            for(int i = 0; i < auranum; i++){
                auradata[i] = DEAD;
            }
            Fade = 1;
            for(int mn = 0; mn < m1num; mn++){
                m1[mn].respawn();
            }
            bulletnum = 0;
        }
        
        
        
        // ------------------------------------
        
        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        
        /* Poll for and process events */
        glfwPollEvents();
    }
    
    glfwTerminate();
    return 0;
}
