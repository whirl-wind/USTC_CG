//88-Line 2D Moving Least Squares Material Point Method (MLS-MPM)[with comments]
#define TC_IMAGE_IO   // Uncomment this line for image exporting functionality
#include <taichi.h>    // Note: You DO NOT have to install taichi or taichi_mpm.
using namespace taichi;// You only need [taichi.h] - see below for instructions.
const int n = 80 /*grid resolution (cells)*/, window_size = 800;
const real dt = 1e-5_f, frame_dt = 1e-3_f, dx = 1.0_f / n, inv_dx = 1.0_f / dx;
auto particle_mass = 1.0_f, vol = 1.0_f;
auto hardening = 10.0_f, E = 1e4_f, nu = 0.2_f;
real mu_0 = E / (2 * (1 + nu)), lambda_0 = E * nu / ((1 + nu) * (1 - 2 * nu));
using Vec = Vector2; using Mat = Matrix2; bool plastic = true;
float vel_strength = 8.0;
bool move_t = false;
/***********************************(1)*****************************************/
struct Particle {
    Vec x, v; Mat F, C; real Jp; int c/*color*/; real mass;
    int ptype/*0: fluid 1: jelly 2: snow*/;
    Particle(Vec x, int c, Vec v = Vec(0), int ptype = 2) : x(x), v(v), F(1), C(0), Jp(1), c(c), ptype(ptype), mass(particle_mass) {}
};
////////////////////////////////////////////////////////////////////////////////
std::vector<Particle> particles;
Vector3 grid[n + 1][n + 1];          // velocity + mass, node_res = cell_res + 1


Vec move_object = Vec(0.96, 0);
real indicate_m = 1.0;
void move_board() {
    if (move_object[0] >= 0.95) {
        indicate_m *= -1.0;
        move_object[1] = 0.0;
    }
    else if (move_object[0] <= 0.05) {
        indicate_m *= -1.0;
        move_object[1] = 0.0;
    }
    if (abs(move_object[1]) <= 500)
        move_object[1] += indicate_m * 5.0;
    else{
        indicate_m *= -1.0;
        move_object[1] = 0.0;
    }

}


bool trans = false;
real indicate = 1.0;
real a = -0.975_f;
int indicate_c = 0xFFFFFF;
void advance(real dt) {
    trans = false;
    indicate_c = 0xFFFFFF;
    std::memset(grid, 0, sizeof(grid));                              // Reset grid
    for (auto& p : particles) {                                             // P2G
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5_f)).cast<int>();//element-wise floor
        Vec fx = p.x * inv_dx - base_coord.cast<real>();
        // Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        Vec w[3]{ Vec(0.5) * sqr(Vec(1.5) - fx), Vec(0.75) - sqr(fx - Vec(1.0)),
                 Vec(0.5) * sqr(fx - Vec(0.5)) };
        /***********************************(2)*****************************************/
        auto e = std::exp(hardening * (1.0_f - p.Jp));
        if (p.ptype == 3) e = 20.0;
        if (p.ptype == 1) e = 0.3;
        auto mu = mu_0 * e, lambda = lambda_0 * e;
        if (p.ptype == 0) mu = 0;
        ////////////////////////////////////////////////////////////////////////////////
        real J = determinant(p.F);         //                         Current volume
        Mat r, s; polar_decomp(p.F, r, s); //Polar decomp. for fixed corotated model
        auto stress =                           // Cauchy stress times dt and inv_dx
            -4 * inv_dx * inv_dx * dt * vol * (2 * mu * (p.F - r) * transposed(p.F) + lambda * (J - 1) * J);
        auto affine = stress + p.mass * p.C;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { // Scatter to grid
            auto dpos = (Vec(i, j) - fx) * dx;
            Vector3 mv(p.v * p.mass, p.mass); //translational momentum
            grid[base_coord.x + i][base_coord.y + j] +=
                w[i].x * w[j].y * (mv + Vector3(affine * dpos, 0));
        }
    }
    for (int i = 0; i <= n; i++) for (int j = 0; j <= n; j++) { //For all grid nodes
        auto& g = grid[i][j];
        if (g[2] > 0) {                                // No need for epsilon here
            g /= g[2];                                   //        Normalize by mass
            g += dt * Vector3(0, -200, 0);               //                  Gravity
            real boundary = 0.05, x = (real)i / n, y = real(j) / n; //boundary thick.,node coord
            if (x < boundary || x > 1 - boundary || y > 1 - boundary) g = Vector3(0); //Sticky
            if (y < boundary) g[1] = std::max(0.0_f, g[1]);             //"Separate"
        }
    }
    for (auto& p : particles) {                                // Grid to particle
        Vector2i base_coord = (p.x * inv_dx - Vec(0.5_f)).cast<int>();//element-wise floor
        Vec fx = p.x * inv_dx - base_coord.cast<real>();
        Vec w[3]{ Vec(0.5) * sqr(Vec(1.5) - fx), Vec(0.75) - sqr(fx - Vec(1.0)),
                 Vec(0.5) * sqr(fx - Vec(0.5)) };
        p.C = Mat(0); p.v = Vec(0);
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
            auto dpos = (Vec(i, j) - fx),
                grid_v = Vec(grid[base_coord.x + i][base_coord.y + j]);
            auto weight = w[i].x * w[j].y;
            p.v += weight * grid_v;                                      // Velocity
            if (p.ptype == 3&& move_t) {
                p.v += Vec(move_object[1], 0.0);
                p.v[1] = 0.0;
            }
            p.C += 4 * inv_dx * Mat::outer_product(weight * grid_v, dpos); // APIC C
        }
        p.x += dt * p.v;                                                // Advection
        auto F = (Mat(1) + dt * p.C) * p.F;                      // MLS-MPM F-update
    /***********************************(3)*****************************************/
        if (p.ptype == 0) { p.F = Mat(1) * sqrt(determinant(F)); }
        else if (p.ptype == 1) { p.F = F; }
        else if (p.ptype == 2) {
            Mat svd_u, sig, svd_v; svd(F, svd_u, sig, svd_v);
            for (int i = 0; i < 2 * int(plastic); i++)                // Snow Plasticity
                sig[i][i] = clamp(sig[i][i], 1.0_f - 2.5e-2_f, 1.0_f + 7.5e-3_f);
            real oldJ = determinant(F); F = svd_u * sig * transposed(svd_v);
            real Jp_new = clamp(p.Jp * oldJ / determinant(F), 0.6_f, 20.0_f);
            p.Jp = Jp_new; p.F = F;
        }
        else if (p.ptype == 3) {
            //std::cout << indicate << std::endl;
            /*Mat svd_u, sig, svd_v; svd(F, svd_u, sig, svd_v);
            for (int i = 0; i < 2 * int(plastic); i++)
                sig[i][i] = clamp(sig[i][i] * (0.0_f - indicate * 5), 0.9729_f , 1.025_f);
            real oldJ = determinant(F); F = svd_u * sig * transposed(svd_v);
            real Jp_new = clamp(p.Jp * oldJ / determinant(F), 0.6_f, 20.0_f);
            p.Jp = Jp_new; p.F = F;*/
            p.F = F;
        }
        ////////////////////////////////////////////////////////////////////////////////
    }
    
    if (!trans) {
        if (indicate > 0) {
            indicate *= a;
        }
        else {
            indicate *= 1 / a;
        }
    }
    
}
/***********************************(4)*****************************************/
void add_object(Vec center, int c, real m= particle_mass, int ptype = 2) {   // Seed particles with position and color
    for (int i = 0; i < 1000; i++)  // Randomly sample 1000 particles in the square
    {
        auto p = Particle((Vec::rand() * 2.0_f - Vec(1)) * 0.08_f + center, c, Vec(0.0, 0.0), ptype);
        p.mass = m;
        particles.push_back(p);
    }
}
////////////////////////////////////////////////////////////////////////////////
/***********************************(5)*****************************************/
void add_object_triangle(Vec v1, Vec v2, Vec v3, int c, int num = 500, Vec velocity = Vec(0.0_f), real m = particle_mass, int ptype = 2)
{
    Vec a1 = v2 - v3;
    Vec a2 = v3 - v1;
    Vec a3 = v1 - v2;

    int i = 0;
    while (i < num) {
        auto pos = Vec::rand();
        real i1 = cross(a2, pos - v1);
        real i2 = cross(a3, pos - v2);
        real i3 = cross(a1, pos - v3);
        if (i1*i2>0 && i2*i3>0&& i3*i1>0) {
            auto p =Particle(pos, c, velocity, ptype);
            p.mass = m;
            particles.push_back(p);
            i++;
        }
    }
}
void add_object_rectangle(Vec v1, Vec v2, int c, int num = 500, Vec velocity = Vec(0.0_f), real m = particle_mass, int ptype = 2)
{
    Vec box_min(min(v1.x, v2.x), min(v1.y, v2.y)), box_max(max(v1.x, v2.x), max(v1.y, v2.y));
    int i = 0;
    while (i < num) {
        auto pos = Vec::rand();
        if (pos.x > box_min.x && pos.y > box_min.y && pos.x < box_max.x && pos.y < box_max.y) {
            auto p = Particle(pos, c, velocity, ptype);
            p.mass = m;
            particles.push_back(p);
            i++;
        }
    }
}
void add_jet(int ptype = 2) {
    add_object_rectangle(Vec(0.05, 0.5), Vec(0.06, 0.51), 0x99FF99, 10, Vec(7.0, 0.0), ptype);
    //add_object_rectangle(Vec(0.5, 0.5), Vec(0.51, 0.51), 0x87CEFA, 10, Vec(0.0, -10.0));
}
void fill_water(float h) {
    add_object_rectangle(Vec(0.05, 0.05), Vec(0.95/*0.78*/, h), 0x87CEFA, h*10000, Vec(0.0, 0.0), 1, 0);
}
void interact() {

}

////////////////////////////////////////////////////////////////////////////////
int main() {
    GUI gui("Real-time 2D MLS-MPM", window_size, window_size);

    /***********************************(6)*****************************************/
    //add_object(Vec(0.86, 0.61), 0x87CEFA, 1, 0);
    add_object(Vec(0.70, 0.80), 0xFFFFCC, 0.9, 2);
    //add_object(Vec(0.86, 0.29), 0x87CEFA, 1, 0);
    //add_object(Vec(0.20, 0.70), 0xED553B, 1, 1);
    add_object_rectangle(Vec(0.12, 0.70), Vec(0.17, 0.75), 0xFF6666, 150, Vec(0.0, 0.0), 0.1, 1);
    add_object_rectangle(Vec(0.22, 0.70), Vec(0.27, 0.75), 0xFF3333, 150, Vec(0.0, 0.0), 1.0, 1);
    add_object_rectangle(Vec(0.32, 0.70), Vec(0.37, 0.75), 0xCC0000, 150, Vec(0.0, 0.0), 10.0,1);
    fill_water(0.4);
    for (int i = 0; i < 3; i++) {
        add_object_triangle(Vec(0.23 * (i + 1) - 0.03, 0.04), Vec(0.23 * (i+1), 0.1), Vec(0.23 * (i + 1) + 0.03, 0.04), 0xFF6600, 500, Vec(0.0_f, 30), 100, 2);
    }
    //add_object_triangle(Vec(0.80, 0.04), Vec(0.88, 0.1), Vec(0.96, 0.04), 0xFF6600, 500, Vec(0.0_f),100, 3);
    //add_object(Vec(0.87, 0.12), 0xED553B, 3);
    add_object_rectangle(Vec(0.05, 0.60), Vec(0.95, 0.65), 0x99FF99, 1000, Vec(0.0, 5.0), 1, 1);
    ////////////////////////////////////////////////////////////////////////////////
    auto& canvas = gui.get_canvas(); int f = 0;
    for (int i = 0;; i++) {                              //              Main Loop
        move_board();
        advance(dt);                                       //     Advance simulation
        if (i % int(frame_dt / dt) == 0) {                 //        Visualize frame
            canvas.clear(0x112F41);                          //       Clear background
            canvas.rect(Vec(0.04), Vec(0.96)).radius(2).color(0x4FB99F).close();// Box
            for (auto p : particles)canvas.circle(p.x).radius(2).color(p.c);//Particles
            gui.update();                                              // Update image
            canvas.img.write_as_image(fmt::format("tmp/{:05d}.png", f++));
            /***********************************(7)*****************************************/
            //if (i < 5e3) add_jet(0);
            ////////////////////////////////////////////////////////////////////////////////
        }
    }
} //----------------------------------------------------------------------------
//This sample shows how to include simulation of water(fluid) & jello(elastic 
//object). Follow the comment marks in the code to see where to make changes
//
//(1) Change particle struct to allow type selection
//(2)&(3) Adjust F update schemes and mpm parameters for different materials
//(4) Adjust add_object to allow type selection
//(5) Sample particles with different materials 