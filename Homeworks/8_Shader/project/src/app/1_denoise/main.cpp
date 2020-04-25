#include <UGL/UGL>
#include <UGM/UGM>

#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "../../tool/Camera.h"
#include "../../tool/SimpleLoader.h"

#include <iostream>

using namespace Ubpa;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow *window);
gl::Texture2D loadTexture(char const* path);
gl::Texture2D genDisplacementmap(const SimpleLoader::OGLResources* resources);

// settings
unsigned int scr_width = 800;
unsigned int scr_height = 600;
float displacement_bias = 0.f;
float displacement_scale = 1.f;
float displacement_lambda = 0.2f;
bool have_denoise = false;

// camera
Camera camera(pointf3(0.0f, 0.0f, 3.0f));
float lastX = scr_width / 2.0f;
float lastY = scr_height / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(scr_width, scr_height, "HW8 - denoise", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    gl::Enable(gl::Capability::DepthTest);

    // build and compile our shader zprogram
    // ------------------------------------
    gl::Shader vs(gl::ShaderType::VertexShader, "../data/shaders/p3t2n3_denoise.vert");
    gl::Shader fs(gl::ShaderType::FragmentShader, "../data/shaders/light.frag");
    gl::Program program(&vs, &fs);
    rgbf ambient{ 0.2f,0.2f,0.2f };
    program.SetTex("albedo_texture", 0);
    program.SetTex("displacementmap", 1);
    program.SetVecf3("point_light_pos", { 0,5,0 });
    program.SetVecf3("point_light_radiance", { 100,100,100 });
    program.SetVecf3("ambient_irradiance", ambient);
    program.SetFloat("roughness", 0.5f );
    program.SetFloat("metalness", 0.f);

    // load model
    // ------------------------------------------------------------------
    auto spot = SimpleLoader::LoadObj("../data/models/spot_triangulated_good.obj", true);
    // world space positions of our cubes
    pointf3 instancePositions[] = {
        pointf3(0.0f,  0.0f,  0.0f),
        pointf3(2.0f,  5.0f, -15.0f),
        pointf3(-1.5f, -2.2f, -2.5f),
        pointf3(-3.8f, -2.0f, -12.3f),
        pointf3(2.4f, -0.4f, -3.5f),
        pointf3(-1.7f,  3.0f, -7.5f),
        pointf3(1.3f, -2.0f, -2.5f),
        pointf3(1.5f,  2.0f, -2.5f),
        pointf3(1.5f,  0.2f, -1.5f),
        pointf3(-1.3f,  1.0f, -1.5f)
    };

    // load and create a texture 
    // -------------------------
    gl::Texture2D spot_albedo = loadTexture("../data/textures/spot_albedo.png");

    gl::Texture2D displacementmap = genDisplacementmap(spot);

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        gl::ClearColor({ ambient, 1.0f });
        gl::Clear(gl::BufferSelectBit::ColorBufferBit | gl::BufferSelectBit::DepthBufferBit); // also clear the depth buffer now!

        program.SetVecf3("camera_pos", camera.Position);

        // bind textures on corresponding texture units
        program.Active(0, &spot_albedo);
        program.Active(1, &displacementmap);

        // pass projection matrix to shader (note that in this case it could change every frame)
        transformf projection = transformf::perspective(to_radian(camera.Zoom), (float)scr_width / (float)scr_height, 0.1f, 100.f);
        program.SetMatf4("projection", projection);

        // camera/view transformation
        program.SetMatf4("view", camera.GetViewMatrix());
        program.SetFloat("displacement_bias", displacement_bias);
        program.SetFloat("displacement_scale", displacement_scale);
        program.SetFloat("displacement_lambda", displacement_lambda);
        program.SetBool("have_denoise", have_denoise);

        // render spots
        for (unsigned int i = 0; i < 10; i++)
        {
            // calculate the model matrix for each object and pass it to shader before drawing
            float angle = 20.0f * i + 10.f * (float)glfwGetTime();
            transformf model(instancePositions[i], quatf{ vecf3(1.0f, 0.3f, 0.5f), to_radian(angle) });
            program.SetMatf4("model", model);
            spot->va->Draw(&program);
        }

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    delete spot;

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(Camera::Movement::FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(Camera::Movement::BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(Camera::Movement::LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(Camera::Movement::RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        camera.ProcessKeyboard(Camera::Movement::UP, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        camera.ProcessKeyboard(Camera::Movement::DOWN, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        have_denoise = !have_denoise;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    gl::Viewport({ 0, 0 }, width, height);
    scr_width = width;
    scr_height = height;
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = static_cast<float>(xpos);
        lastY = static_cast<float>(ypos);
        firstMouse = false;
    }

    float xoffset = static_cast<float>(xpos) - lastX;
    float yoffset = lastY - static_cast<float>(ypos); // reversed since y-coordinates go from bottom to top

    lastX = static_cast<float>(xpos);
    lastY = static_cast<float>(ypos);

    camera.ProcessMouseMovement(static_cast<float>(xoffset), static_cast<float>(yoffset));
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}

gl::Texture2D loadTexture(char const* path)
{
    gl::Texture2D tex;
    tex.SetWrapFilter(gl::WrapMode::Repeat, gl::WrapMode::Repeat, gl::MinFilter::Linear, gl::MagFilter::Linear);
    // load image, create texture and generate mipmaps
    int width, height, nrChannels;
    stbi_set_flip_vertically_on_load(true); // tell stb_image.h to flip loaded texture's on the y-axis.
    unsigned char* data = stbi_load(path, &width, &height, &nrChannels, 0);
    gl::PixelDataFormat c2f[4] = {
        gl::PixelDataFormat::Red,
        gl::PixelDataFormat::Rg,
        gl::PixelDataFormat::Rgb,
        gl::PixelDataFormat::Rgba
    };
    gl::PixelDataInternalFormat c2if[4] = {
        gl::PixelDataInternalFormat::Red,
        gl::PixelDataInternalFormat::Rg,
        gl::PixelDataInternalFormat::Rgb,
        gl::PixelDataInternalFormat::Rgba
    };
    if (data)
    {
        tex.SetImage(0, c2if[nrChannels - 1], width, height, c2f[nrChannels - 1], gl::PixelDataType::UnsignedByte, data);
        tex.GenerateMipmap();
    }
    else
    {
        std::cout << "Failed to load texture" << std::endl;
    }
    stbi_image_free(data);

    return tex;
}

gl::Texture2D genDisplacementmap(const SimpleLoader::OGLResources* resources) {
    float* displacementData = new float[1024 * 1024];
    // TODO: HW8 - 1_denoise | genDisplacementmap
    // 1. set displacementData with resources's positions, indices, normals, ...
    size_t vertexNum = resources->positions.size();
    std::vector<pointf3> new_positions;
    std::vector<int> neighbour_num;
    std::vector<float> delta;
    bool* neighbour = new bool[vertexNum * vertexNum];
    new_positions.resize(vertexNum);
    neighbour_num.resize(vertexNum);
    delta.resize(vertexNum);
    for (size_t i = 0; i < vertexNum; i++) {
        new_positions[i] = pointf3(0.0, 0.0, 0.0);
        neighbour_num[i] = 0;
        for (size_t j = 0; j < vertexNum; j++) {
            neighbour[j * (vertexNum - 1) + i] = false;
        }
    }
    for (size_t i = 0; i < resources->indices.size(); i += 3) {
        auto i1 = resources->indices[i];
        auto i2 = resources->indices[i + 1];
        auto i3 = resources->indices[i + 2];

        pointf3 v1 = resources->positions.at(i1);
        pointf3 v2 = resources->positions.at(i2);
        pointf3 v3 = resources->positions.at(i3);

        if (!neighbour[i1 * (vertexNum - 1) + i2]) {
            for (size_t j = 0; j < 3; j++) {
                new_positions[i1].at(j) += v2.at(j);
                new_positions[i2].at(j) += v1.at(j);
            }
            neighbour[i1 * (vertexNum - 1) + i2] = true;
            neighbour[i2 * (vertexNum - 1) + i1] = true;
            neighbour_num[i1] += 1;
            neighbour_num[i2] += 1;
        }
        if (!neighbour[i2 * (vertexNum - 1) + i3]) {
            for (size_t j = 0; j < 3; j++) {
                new_positions[i2].at(j) += v3.at(j);
                new_positions[i3].at(j) += v2.at(j);
            }
            neighbour[i2 * (vertexNum - 1) + i3] = true;
            neighbour[i3 * (vertexNum - 1) + i2] = true;
            neighbour_num[i2] += 1;
            neighbour_num[i3] += 1;
        }
        if (!neighbour[i3 * (vertexNum - 1) + i1]) {
            for (size_t j = 0; j < 3; j++) {
                new_positions[i3].at(j) += v1.at(j);
                new_positions[i1].at(j) += v3.at(j);
            }
            neighbour[i3 * (vertexNum - 1) + i1] = true;
            neighbour[i1 * (vertexNum - 1) + i3] = true;
            neighbour_num[i3] += 1;
            neighbour_num[i1] += 1;
        }
    }

    vecf3 d0;
    for (size_t j = 0; j < 3; j++) {
        d0.at(j) = resources->positions[0].at(j) - new_positions[0].at(j) / float(neighbour_num[0]);
    }
    delta[0] = d0.dot(resources->normals.at(0).cast_to<vecf3>());
    float max = delta[0], min = delta[0];
    for (size_t i = 1; i < vertexNum; i++) {
        vecf3 d;
        for (size_t j = 0; j < 3; j++) {
            d.at(j) = resources->positions[i].at(j) - new_positions[i].at(j) / float(neighbour_num[i]);
        }
        delta[i] = d.dot(resources->normals.at(i).cast_to<vecf3>());
        if (delta[i] < min) min = delta[i];
        else if (delta[i] > max) max = delta[i];
    }

    // 2. change global variable: displacement_bias, displacement_scale, displacement_lambda
    displacement_scale = (max - min);
    displacement_lambda = 0.75f;
    displacement_bias = min;
    //std::cout << min << " " << max << std::endl;
    for (size_t i = 0; i < resources->indices.size(); i += 3) {
        auto i1 = resources->indices[i];
        auto i2 = resources->indices[i + 1];
        auto i3 = resources->indices[i + 2];
        
        int u1 = (int)std::round(1024 * std::clamp(resources->texcoords.at(i1).at(0), 0.f, 1.f) - 0.5);
        int v1 = (int)std::round(1024 * std::clamp(resources->texcoords.at(i1).at(1), 0.f, 1.f) - 0.5);
        int u2 = (int)std::round(1024 * std::clamp(resources->texcoords.at(i2).at(0), 0.f, 1.f) - 0.5);
        int v2 = (int)std::round(1024 * std::clamp(resources->texcoords.at(i2).at(1), 0.f, 1.f) - 0.5);
        int u3 = (int)std::round(1024 * std::clamp(resources->texcoords.at(i3).at(0), 0.f, 1.f) - 0.5);
        int v3 = (int)std::round(1024 * std::clamp(resources->texcoords.at(i3).at(1), 0.f, 1.f) - 0.5);
        //v1 = 1023 - v1;
        //v2 = 1023 - v2;
        //v3 = 1023 - v3;

        float a1 = delta[i1];
        float a2 = delta[i2];
        float a3 = delta[i3];

        displacementData[u1 + 1024 * v1] = (a1 - displacement_bias) / displacement_scale;
        displacementData[u2 + 1024 * v2] = (a2 - displacement_bias) / displacement_scale;
        displacementData[u3 + 1024 * v3] = (a3 - displacement_bias) / displacement_scale;
        
        int bu_min = u1, bu_max = u1, bv_min = v1, bv_max = v1;
        if (bu_min > u2) bu_min = u2;
        else if (bu_max < u2) bu_max = u2;
        if (bu_min > u3) bu_min = u3;
        else if (bu_max < u3) bu_max = u3;
        if (bv_min > v2) bv_min = v2;
        else if (bv_max < v2) bv_max = v2;
        if (bv_min > v3) bv_min = v3;
        else if (bv_max < v3) bv_max = v3;

        for (int u = bu_min; u <= bu_max; u++) {
            for (int v = bv_min; v <= bv_max; v++) {
                float d3 = (u - u1) * (v2 - v1) - (v - v1) * (u2 - u1);
                float d1 = (u - u2) * (v3 - v2) - (v - v2) * (u3 - u2);
                float d2 = (u - u3) * (v1 - v3) - (v - v3) * (u1 - u3);
                if (d1 * d2 >= 0.0 && d2 * d3 >= 0.0) {
                    float d0 = d1 + d2 + d3;
                    displacementData[u + 1024 * v] = (a1 * d1 / d0 + a2 * d2 / d0 + a3 * d3 / d0 - displacement_bias) / displacement_scale;
                }
            }
        }
    }

    gl::Texture2D displacementmap;
    displacementmap.SetImage(0, gl::PixelDataInternalFormat::Red, 1024, 1024, gl::PixelDataFormat::Red, gl::PixelDataType::Float, displacementData);
    displacementmap.SetWrapFilter(gl::WrapMode::Repeat, gl::WrapMode::Repeat,
        gl::MinFilter::Linear, gl::MagFilter::Linear);
    stbi_uc* stbi_data = new stbi_uc[1024 * 1024];
    for (size_t i = 0; i < 1024 * 1024; i++)
        stbi_data[i] = static_cast<stbi_uc>(std::clamp(displacementData[i] * 255.f, 0.f, 255.f));
    stbi_write_png("../data/1_denoise_displacement_map.png", 1024, 1024, 1, stbi_data, 1024);
    delete[] stbi_data;
    delete[] displacementData;
    return displacementmap;
}
