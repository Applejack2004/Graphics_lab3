#version 430

out vec4 FragColor;
in vec3 glPosition;

/*** DATA STRUCTURES ***/
struct SCamera {
    vec3 Position;
    vec3 View;
    vec3 Up;
    vec3 Side;
    // aspect ratio
    vec2 Scale;
};

struct SRay {
    vec3 Origin;
    vec3 Direction;
    float RefractIndice;
};

struct STracingRay
{
    SRay Ray;
    float Contrib;
    int Depth;
};

struct SIntersection {
    float Time;
    vec3 Point;
    vec3 Normal;
    vec3 Color;
    // ambient, diffuse and specular coeffs
    vec4 LightCoeffs;
    float ReflectionCoef;
    float RefractionCoef;
    float RefractIndice;
};

struct SMaterial {
    vec3 Color;
    // ambient, diffuse and specular coeffs
    vec4 LightCoeffs;
    float ReflectionCoef;
    float RefractionCoef;
    float RefractIndice;
};

struct SLight
{
    vec3 Position;
    vec3 Color;
    float Intensity;
};

struct STriangle {
    vec3 v1;
    vec3 v2;
    vec3 v3;
    int MaterialIdx;
};

struct SSphere {
    vec3 Center;
    float Radius;
    int MaterialIdx;
};


/*** RAY STACK ***/
STracingRay Stack[100];
int StackSize = 0;

void PushRay(STracingRay ray) {
    Stack[StackSize++] = ray;
}

STracingRay PopRay() {
    return Stack[--StackSize];
}

bool IsEmpty() {
    return StackSize <= 0;
}

/*** GLOBALS ***/
#define EPSILON 0.001
#define BIG 1000000.0

int TraceDepth = 100;

SCamera uCamera;
SMaterial materials[8];
SLight lights[2];
STriangle triangles[16];
SSphere spheres[2];

/*** INTERSECTIONS ***/
bool IntersectSphere(SSphere sphere, SRay ray, float start, float final, out float
    time)
{
    ray.Origin -= sphere.Center;
    float A = dot(ray.Direction, ray.Direction);
    float B = dot(ray.Direction, ray.Origin);
    float C = dot(ray.Origin, ray.Origin) - sphere.Radius * sphere.Radius;
    float D = B * B - A * C;
    if (D > 0.0)
    {
        D = sqrt(D);
        //time = min ( max ( 0.0, ( -B - D ) / A ), ( -B + D ) / A );
        float t1 = (-B - D) / A;
        float t2 = (-B + D) / A;
        if (t1 < 0 && t2 < 0)
            return false;

        if (min(t1, t2) < 0)
        {
            time = max(t1, t2);
            return true;
        }
        time = min(t1, t2);
        return true;
    }
    return false;
}


bool IntersectTriangle(SRay ray, vec3 v1, vec3 v2, vec3 v3, out float time)
{
    // // Compute the intersection of ray with a triangle using geometric solution
    // Input: // points v0, v1, v2 are the triangle's vertices
    // rayOrig and rayDir are the ray's origin (point) and the ray's direction
    // Return: // return true is the ray intersects the triangle, false otherwise
    // bool intersectTriangle(point v0, point v1, point v2, point rayOrig, vector rayDir) {
    // compute plane's normal vector
    time = -1;
    vec3 A = v2 - v1;
    vec3 B = v3 - v1;
    // no need to normalize vector
    vec3 N = cross(A, B);
    // N
    // // Step 1: finding P
    44
        // // check if ray and plane are parallel ?
        float NdotRayDirection = dot(N, ray.Direction);
    if (abs(NdotRayDirection) < 0.001)
        return false;
    // they are parallel so they don't intersect !
    // compute d parameter using equation 2
    float d = dot(N, v1);
    // compute t (equation 3)
    float t = -(dot(N, ray.Origin) - d) / NdotRayDirection;
    // check if the triangle is in behind the ray
    if (t < 0)
        return false;
    // the triangle is behind
    // compute the intersection point using equation 1
    vec3 P = ray.Origin + t * ray.Direction;
    // // Step 2: inside-outside test //
    vec3 C;
    // vector perpendicular to triangle's plane
    // edge 0
    vec3 edge1 = v2 - v1;
    vec3 VP1 = P - v1;
    C = cross(edge1, VP1);
    if (dot(N, C) < 0)
        return false;
    // P is on the right side
    // edge 1
    vec3 edge2 = v3 - v2;
    vec3 VP2 = P - v2;
    C = cross(edge2, VP2);
    if (dot(N, C) < 0)
        return false;
    // P is on the right side
    // edge 2
    vec3 edge3 = v1 - v3;
    vec3 VP3 = P - v3;
    C = cross(edge3, VP3);
    if (dot(N, C) < 0)
        return false;
    // P is on the right side;
    time = t;
    return true;
    // this ray hits the triangle
}

/*** FUNCTIONS ***/
SCamera InitializeDefaultCamera() {
    //** CAMERA **//
    SCamera camera;
    camera.Position = vec3(1.0, 0.0, -4.99);
    camera.View = vec3(0.0, 0.0, 1.0);
    camera.Up = vec3(0.0, 1.0, 0.0);
    camera.Side = vec3(1.778, 0.0, 0.0);
    float mult = 1.0;
    camera.Scale = vec2(mult * 1.0, mult * 1.0);
    return camera;
}

void InitializeDefaultScene() {

    /** TRIANGLES **/
    /* left wall */
    triangles[0].v1 = vec3(-5.0, -5.0, -5.0);
    triangles[0].v2 = vec3(-5.0, 5.0, 5.0);
    triangles[0].v3 = vec3(-5.0, 5.0, -5.0);
    triangles[0].MaterialIdx = 0;  
    triangles[1].v1 = vec3(-5.0, -5.0, -5.0);
    triangles[1].v2 = vec3(-5.0, -5.0, 5.0);
    triangles[1].v3 = vec3(-5.0, 5.0, 5.0);
    triangles[1].MaterialIdx = 0; 

    /* right wall */
    triangles[2].v1 = vec3(5.0, -5.0, -5.0);
    triangles[2].v2 = vec3(5.0, 5.0, -5.0);
    triangles[2].v3 = vec3(5.0, -5.0, 5.0);
    triangles[2].MaterialIdx = 1; 
    triangles[3].v1 = vec3(5.0, 5.0, 5.0);
    triangles[3].v2 = vec3(5.0, -5.0, 5.0);
    triangles[3].v3 = vec3(5.0, 5.0, -5.0);
    triangles[3].MaterialIdx = 1;  

    /* up wall */
    triangles[4].v1 = vec3(-5.0, 5.0, -5.0);
    triangles[4].v2 = vec3(-5.0, 5.0, 5.0);
    triangles[4].v3 = vec3(5.0, 5.0, -5.0);
    triangles[4].MaterialIdx = 2; 
    triangles[5].v1 = vec3(5.0, 5.0, 5.0);
    triangles[5].v2 = vec3(5.0, 5.0, -5.0);
    triangles[5].v3 = vec3(-5.0, 5.0, 5.0);
    triangles[5].MaterialIdx = 2;  

    /* down wall */
    triangles[6].v1 = vec3(-5.0, -5.0, -5.0);
    triangles[6].v2 = vec3(5.0, -5.0, -5.0);
    triangles[6].v3 = vec3(-5.0, -5.0, 5.0);
    triangles[6].MaterialIdx = 3;  
    triangles[7].v1 = vec3(5.0, -5.0, 5.0);
    triangles[7].v2 = vec3(-5.0, -5.0, 5.0);
    triangles[7].v3 = vec3(5.0, -5.0, -5.0);
    triangles[7].MaterialIdx = 3;  

    /* front wall */
    triangles[8].v1 = vec3(-5.0, -5.0, 5.0);
    triangles[8].v2 = vec3(5.0, -5.0, 5.0);
    triangles[8].v3 = vec3(-5.0, 5.0, 5.0);
    triangles[8].MaterialIdx = 4; 
    triangles[9].v1 = vec3(5.0, 5.0, 5.0);
    triangles[9].v2 = vec3(-5.0, 5.0, 5.0);
    triangles[9].v3 = vec3(5.0, -5.0, 5.0);
    triangles[9].MaterialIdx = 4; 

    /* back wall */
    triangles[10].v1 = vec3(-5.0, -5.0, -5.0);
    triangles[10].v2 = vec3(5.0, -5.0, -5.0);
    triangles[10].v3 = vec3(-5.0, 5.0, -5.0);
    triangles[10].MaterialIdx = 5;  
    triangles[11].v1 = vec3(5.0, 5.0, -5.0);
    triangles[11].v2 = vec3(-5.0, 5.0, -5.0);
    triangles[11].v3 = vec3(5.0, -5.0, -5.0);
    triangles[11].MaterialIdx = 5;  

    /** SPHERES **/
   

    /** TETRAHEDRON **/
   

void InitializeDefaultLightMaterials() {
    //** LIGHTS **//
    lights[0].Position = vec3(4.0, 3.0, -2.0);
    lights[0].Color = vec3(1.0, 1.0, 0.82);
    lights[0].Intensity = 0.6;

    lights[1].Position = vec3(-3.0, -2.0, -4.0);
    lights[1].Color = vec3(1.0, 0.85, 0.7);
    lights[1].Intensity = 0.3;

    /** MATERIALS **/
    vec4 lightCoeffs = vec4(0.4, 0.9, 0.2, 512.0);
    materials[0].Color = vec3(0.715, 0.109, 0.109);  // Red
    materials[0].LightCoeffs = lightCoeffs;
    materials[0].ReflectionCoef = 0.0;
    materials[0].RefractionCoef = 0.0;
    materials[0].RefractIndice = 1.0;

    materials[1].Color = vec3(0.645, 0.836, 0.652);  // Green
    materials[1].LightCoeffs = lightCoeffs;
    materials[1].ReflectionCoef = 0.0;
    materials[1].RefractionCoef = 0.0;
    materials[1].RefractIndice = 1.0;

    materials[2].Color = vec3(0.621, 0.656, 0.852);  // Blue
    materials[2].LightCoeffs = lightCoeffs;
    materials[2].ReflectionCoef = 0.0;
    materials[2].RefractionCoef = 0.0;
    materials[2].RefractIndice = 1.0;

    materials[3].Color = vec3(1.0, 1.0, 0.0);  // Yellow
    materials[3].LightCoeffs = lightCoeffs;
    materials[3].ReflectionCoef = 0.0;
    materials[3].RefractionCoef = 0.0;
    materials[3].RefractIndice = 1.0;

    materials[4].Color = vec3(0.0, 1.0, 1.0);  // Cyan
    materials[4].LightCoeffs = lightCoeffs;
    materials[4].ReflectionCoef = 0.0;
    materials[4].RefractionCoef = 0.0;
    materials[4].RefractIndice = 1.0;

    materials[5].Color = vec3(1.0, 0.0, 1.0);  // Purple
    materials[5].LightCoeffs = lightCoeffs;
    materials[5].ReflectionCoef = 0.0;
    materials[5].RefractionCoef = 0.0;
    materials[5].RefractIndice = 1.0;

    materials[6].Color = vec3(1.0, 1.0, 1.0);
    materials[6].LightCoeffs = vec4(0.2, 0.5, 0.9, 1024.0);
    materials[6].ReflectionCoef = 0.8;
    materials[6].RefractionCoef = 0.0;
    materials[6].RefractIndice = 1.0;

    materials[7].Color = vec3(1.0, 1.0, 1.0);
    materials[7].LightCoeffs = vec4(0.2, 0.5, 0.9, 16.0);
    materials[7].ReflectionCoef = 0.0;
    materials[7].RefractionCoef = 0.7;
    materials[7].RefractIndice = 2.0;
}

SRay GenerateRay() {
    vec2 coords = glPosition.xy * uCamera.Scale;
    vec3 direction = uCamera.View + uCamera.Side * coords.x + uCamera.Up * coords.y;
    return SRay(uCamera.Position, normalize(direction), 1.0);
}

bool Raytrace(SRay ray, SSphere spheres[], STriangle triangles[], SMaterial
    materials[], float start, float final, inout SIntersection intersect)
{
    bool result = false;
    float test = start;
    intersect.Time = final;
    //calculate intersect with spheres
    for (int i = 0; i < 2; i++)
    {
        SSphere sphere = spheres[i];
        if (IntersectSphere(sphere, ray, start, final, test) && test < intersect.Time)
        {
            intersect.Time = test;
            intersect.Point = ray.Origin + ray.Direction * test;
            intersect.Normal = normalize(intersect.Point - spheres[i].Center);
            intersect.Color = vec3(1, 0, 0);
            intersect.LightCoeffs = vec4(0, 0, 0, 0);
            intersect.ReflectionCoef = 0;
            intersect.RefractionCoef = 0;
            intersect.MaterialType = 0;
            result = true;
        }
    }
    //calculate intersect with triangles
    for (int i = 0; i < 10; i++)
    {
        STriangle triangle = triangles[i];
        if (IntersectTriangle(ray, triangle.v1, triangle.v2, triangle.v3, test)
            && test < intersect.Time)
        {
            intersect.Time = test;
            intersect.Point = ray.Origin + ray.Direction * test;
            intersect.Normal =
                normalize(cross(triangle.v1 - triangle.v2, triangle.v3 - triangle.v2));
            intersect.Color = vec3(1, 0, 0);
            intersect.LightCoeffs = vec4(0, 0, 0, 0);
            intersect.ReflectionCoef = 0;
            intersect.RefractionCoef = 0;
            intersect.MaterialType = 0;
            result = true;
        }
    }
    return result;
}


float Shadow(SIntersection intersect, int lightNumber) {
    // Point is lighted
    float shadowing = 1.0;
    // Vector to the light source
    vec3 direction = normalize(lights[lightNumber].Position - intersect.Point);
    // Distance to the light source
    float distanceLight = distance(lights[lightNumber].Position, intersect.Point);
    // Generation shadow ray for this light source
    SRay shadowRay = SRay(intersect.Point + direction * EPSILON, direction, 1.0);
    // ...test intersection this ray with each scene object
    SIntersection shadowIntersect;
    shadowIntersect.Time = BIG;
    // trace ray from shadow ray begining to light source position
    if (Raytrace(shadowRay, 0, distanceLight, shadowIntersect)) {
        // this light source is invisible in the intercection point
        shadowing = 0.0;
    }
    return shadowing;
}

vec3 Phong(SIntersection intersect) {
    vec3 result = vec3(0.0, 0.0, 0.0);
    for (int i = 0; i < 2; ++i) {
        vec3 light = normalize(lights[i].Position - intersect.Point);
        float diffuse = max(dot(light, intersect.Normal), 0.0);
        vec3 view = normalize(uCamera.Position - intersect.Point);
        vec3 reflected = reflect(-view, intersect.Normal);
        float specular = pow(max(dot(reflected, light), 0.0), intersect.LightCoeffs.w);
        result += lights[i].Intensity * (intersect.LightCoeffs.x * intersect.Color * lights[i].Color +
            intersect.LightCoeffs.y * diffuse * intersect.Color * lights[i].Color * Shadow(intersect, i) +
            intersect.LightCoeffs.z * specular);
    }
    return result;
}

void main(void)
{
    uCamera = InitializeDefaultCamera();
    SRay ray = GenerateRay();

    SIntersection intersect;
    intersect.Time = BIG;
    vec3 resultColor = vec3(0, 0, 0);

    InitializeDefaultScene();
    InitializeDefaultLightMaterials();

    STracingRay curRay = STracingRay(ray, 1, 0);
    PushRay(curRay);
    while (StackSize > 0) {
        curRay = PopRay();
        ray = curRay.Ray;
        SIntersection intersect;
        float start = 0;
        float final = BIG;

        if (Raytrace(ray, start, final, intersect)) {
            float contribution = curRay.Contrib * (1.0 - intersect.ReflectionCoef - intersect.RefractionCoef);
            resultColor += contribution * Phong(intersect);

            if (intersect.ReflectionCoef > 0) {
                vec3 reflectDirection = reflect(ray.Direction, intersect.Normal);
                float contribution = curRay.Contrib * intersect.ReflectionCoef;
                STracingRay reflectRay = STracingRay(
                    SRay(intersect.Point + reflectDirection * EPSILON, reflectDirection, 1.0),
                    contribution, curRay.Depth + 1);
                if (reflectRay.Depth < TraceDepth)
                    PushRay(reflectRay);
            }
            if (intersect.RefractionCoef > 0) {
                vec3 refractionDirection = refract(
                    ray.Direction, intersect.Normal, ray.RefractIndice / intersect.RefractIndice);
                float contribution = curRay.Contrib * intersect.RefractionCoef;
                STracingRay refractRay = STracingRay(
                    SRay(intersect.Point + refractionDirection * EPSILON, refractionDirection,
                        intersect.RefractIndice),
                    contribution, curRay.Depth + 1);
                if (refractRay.Depth < TraceDepth)
                    PushRay(refractRay);
            }
        }
    }
    FragColor = vec4(resultColor, 1.0);
}
