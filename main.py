import math
import pygame

FPS = 60

pygame.init()

c = 0.5

acceleration = 0.002 * c

SPF = 1/FPS

relativity = True
relSwitchPressed = False

canvas = pygame.display.set_mode((1600,1000),pygame.RESIZABLE)
pygame.display.set_caption("Relativity Simulation")
screentf = lambda p: (p[0]+canvas.get_width()/2, p[1]+canvas.get_height()/2)
corpus = [(30,50),(0,20),(-30,50),(0,20),(0,-20),(30,20),(0,-20),(-30,20), (0,-20), (0,20)]


def beta(velocity):
    return velocity/c

def gamma(velocity_x, velocity_y):
    return 1/math.sqrt(1-(velocity_x**2 + velocity_y**2)/c**2)

# The returned linear maps map x, y -> x', y' and x, y -> Delta t (both under the assumption t' = 0)
def lTransform(velocity_x, velocity_y=None):
    if velocity_y is None:
        velocity_y = velocity_x[1]
        velocity_x = velocity_x[0]
    bx = beta(velocity_x)
    by = beta(velocity_y)
    g = gamma(velocity_x, velocity_y)
    h = g**2/(1+g)

    m11, m12, m13 = g,          -g*bx/c,    -g*by/c
    m21, m22, m23 = -g*bx*c,    1+h*bx*bx,    h*by*bx
    m31, m32, m33 = -g*by*c,    h*bx*by,    1+h*by*by

    f1, f2, f3 = 1/m11, -m12/m11, -m13/m11

    return (
        lambda p: (
            m21*(f2*p[0] + f3*p[1]) + m22*p[0] + m23*p[1],
            m31*(f2*p[0] + f3*p[1]) + m32*p[0] + m33*p[1]
        ),
        lambda p: f2 * p[0] + f3 * p[1]
    )

def twoVelLTransform(vel1, vel2):
    bx1 = beta(vel1[0])
    by1 = beta(vel1[1])
    bx2 = beta(vel2[0])
    by2 = beta(vel2[1])
    g1 = gamma(*vel1)
    g2 = gamma(*vel2)
    h1 = g1**2/(1+g1)
    h2 = g2**2/(1+g2)

    a11, a12, a13 = g1, -g1 * bx1 / c, -g1 * by1 / c
    a21, a22, a23 = -g1 * bx1 * c, 1 + h1 * bx1 * bx1, h1 * by1 * bx1
    a31, a32, a33 = -g1 * by1 * c, h1 * bx1 * by1, 1 + h1 * by1 * by1

    b11, b12, b13 = g2, -g2 * bx2 / c, -g2 * by2 / c
    b21, b22, b23 = -g2 * bx2 * c, 1 + h2 * bx2 * bx2, h2 * by2 * bx2
    b31, b32, b33 = -g2 * by2 * c, h2 * bx2 * by2, 1 + h2 * by2 * by2

    m11, m12, m13 = dotProduct((a11,a21,a31),(b11,b12,b13)), dotProduct((a12,a22,a32),(b11,b12,b13)), dotProduct((a13,a23,a33),(b11,b12,b13))
    m21, m22, m23 = dotProduct((a11,a21,a31),(b21,b22,b23)), dotProduct((a12,a22,a32),(b21,b22,b23)), dotProduct((a13,a23,a33),(b21,b22,b23))
    m31, m32, m33 = dotProduct((a11,a21,a31),(b31,b32,b33)), dotProduct((a12,a22,a32),(b31,b32,b33)), dotProduct((a13,a23,a33),(b31,b32,b33))

    f1, f2, f3 = 1/m11, -m12/m11, -m13/m11

    return (
        lambda p: (
            m21*(f2*p[0] + f3*p[1]) + m22*p[0] + m23*p[1],
            m31*(f2*p[0] + f3*p[1]) + m32*p[0] + m33*p[1]
        ),
        lambda p: f2 * p[0] + f3 * p[1]
    )

def posTransform(pos):
    return lambda p: (p[0]+pos[0], p[1]+pos[1])

def dotProduct(p1, p2):
    return sum(map(lambda x, y: x*y, p1, p2))

def vecMul(a, vec):
    return a*vec[0], a*vec[1]

def vecAdd(a, b):
    return a[0]+b[0], a[1]+b[1]

def addVel(v, u):
    v2 = v[0]*v[0] + v[1]*v[1]
    if v2 <= 0:
        return u
    alpha = math.sqrt(1-v2/c**2)
    return vecMul(1/(1+dotProduct(v, u)/c**2),vecAdd(vecMul(alpha,u), vecMul(1 + (1-alpha)*dotProduct(v,u)/v2,v)))

def hueToRgb(p, q, t):
  if t < 0: t += 1
  if t > 1: t -= 1
  if t < 1/6: return p + (q - p) * 6 * t
  if t < 1/2: return q
  if t < 2/3: return p + (q - p) * (2/3 - t) * 6
  return p

def hslToRgb(h, s, l):
  if s == 0 :
    r = g = b = l
  else:
    q = l * (1 + s) if l < 0.5 else l + s - l * s
    p = 2 * l - q
    r = hueToRgb(p, q, h + 1/3)
    g = hueToRgb(p, q, h)
    b = hueToRgb(p, q, h - 1/3)

  return int(r * 255), int(g * 255), int(b * 255)


class InertialSystem:
    def __init__(self, vel, obstacles=None, radio_lamps=None, color = (255,0,0)):
        if radio_lamps is None:
            radio_lamps = []
        if obstacles is None:
            obstacles = []
        self.vel = vel
        self.obstacles = obstacles
        self.radio_lamps = radio_lamps
        self.color = color

        self.pos = 0, 0
        self.t = 0

    def drawAndUpdate(self, canvas, vel, dt_prime):
        postf = posTransform(self.pos)



        if relativity:
            summedVel = addVel(vecMul(-1,self.vel), vel)
            ltf, tLTF = twoVelLTransform(vecMul(-1,self.vel), vel)
            lorentzFactor = gamma(*summedVel)
            dt = dt_prime * lorentzFactor
        else:
            summedVel = vecAdd(vecMul(-1,self.vel), vel)
            dt = dt_prime

        self.pos = self.pos[0] + summedVel[0] * dt, self.pos[1] + summedVel[1] * dt
        self.t += dt


        transformedObstacles = map(lambda o: list(map(lambda p: screentf(ltf(postf(p))), o)), self.obstacles) if relativity else map(lambda o: list(map(lambda p: screentf(postf(p)), o)), self.obstacles)

        for obst in transformedObstacles:
            pygame.draw.polygon(canvas, self.color, obst, 5)

        transformedRadioLamps = map(lambda p: screentf(ltf(postf(p))), self.radio_lamps) if relativity else map(
            lambda p: screentf(postf(p)), self.radio_lamps)
        transformedLampTime = list(
            map(lambda p: self.t - tLTF(postf(p)), self.radio_lamps) if relativity else [self.t] * len(self.radio_lamps))

        for rl, rlt in zip(transformedRadioLamps, transformedLampTime):
            pygame.draw.circle(canvas, hslToRgb(-rlt / 10000 - int(-rlt / 10000) + 1, 0.5, 0.5), rl, 10)

fixedObstacles = [
    [(900,0), (950, -50), (950, -100), (850, -100), (950, -130), (875, -130), (950, -160), (900, -160), (950, -190), (925, -190),(1002, -250),(1075, -190), (1050, -190), (1100, -160), (1050,-160), (1125, -130), (1050, -130),(1150,-100),(1050,-100),(1050,-50),(1100,0)]
] + [[(d,100),(d,200),(100+d,200),(100+d,100)] for d in range(-10000,10000,100)]

fixedRadio_lamps = list(zip(range(-200, -40000, -200),[0]*len(list(range(-200, -40000, -200)))))

rocketDistance = 1000
rocketCount = 1000

inSysts = [
    InertialSystem((0,0), fixedObstacles, fixedRadio_lamps),
    InertialSystem((-0.8*c,0),
                   [((100-d,-200),(200-d,-200),(200-d,-300),(100-d,-300)) for d in range(0,rocketCount*rocketDistance,rocketDistance)] +
                   [((200-d,-200),(200-d,-300),(300-d,-250)) for d in range(0,rocketCount*rocketDistance,rocketDistance)],
                   [(133-d,-250) for d in range(0,rocketCount*rocketDistance,rocketDistance)] +
                   [(167-d,-250) for d in range(0,rocketCount*rocketDistance,rocketDistance)],
                   (0,0,255))
]


font = pygame.font.Font(None, 50)

vel = (0,0)
t = 0

clock = pygame.time.Clock()

end = False
while not end:
    dt = clock.tick(FPS)

    canvas.fill((255,255,255))

    pygame.draw.polygon(canvas, (0, 0, 0), list(map(screentf, corpus)), 5)
    pygame.draw.circle(canvas, (0, 0, 0), screentf((0,-40)), 20, 5)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            end = True

    keys = pygame.key.get_pressed()

    left = keys[pygame.K_LEFT] or keys[pygame.K_a]
    right = keys[pygame.K_RIGHT] or keys[pygame.K_d]
    up = keys[pygame.K_UP] or keys[pygame.K_w]
    down = keys[pygame.K_DOWN] or keys[pygame.K_s]

    if keys[pygame.K_r] and not relSwitchPressed:
        relSwitchPressed = True
        relativity = not relativity
    if not keys[pygame.K_r] and relSwitchPressed:
        relSwitchPressed = False

    dvel = int(left-right), int(up-down)
    lenDVel = math.sqrt(dvel[0]**2 + dvel[1]**2)
    if lenDVel != 0:
        dvel = vecMul(acceleration * dt / lenDVel, dvel)
        if keys[pygame.K_LCTRL]: dvel = vecMul(10, dvel)
        if keys[pygame.K_LSHIFT]: dvel = vecMul(0.1, dvel)

    vel = vecMul(0.9+0.095*(lenDVel!=0), addVel(vel, dvel))

    t += dt

    for inSys in inSysts:
        inSys.drawAndUpdate(canvas, vel, dt)

    text = font.render(f'v = 0.{int(1000*math.sqrt(vel[0]**2 + vel[1]**2)/c)}c     t = {int(inSysts[0].t//1000)}s     t\' = {int(t//1000)}s', True, (0, 200, 200))
    textRect = text.get_rect()
    textRect.x = 5
    textRect.y = 5
    canvas.blit(text, textRect)
    if not relativity:
        noRelText = font.render(f'No relativity', True, (200, 0, 0))
        noRelTextRect = noRelText.get_rect()
        noRelTextRect.x = 5
        noRelTextRect.y = 5 + textRect.bottom
        canvas.blit(noRelText, noRelTextRect)

    pygame.display.update()