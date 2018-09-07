from pymol import *

cmd.select("s", "all")
model = cmd.get_model("s", "1")

total_mass = 0.0
x = 0
y = 0
z = 0

for m in model.atom:
    m_mass = m.get_mass()
    x += m.coord[0] * m_mass
    y += m.coord[1] * m_mass
    z += m.coord[2] * m_mass
    total_mass += m_mass

x /= total_mass
y /= total_mass
z /= total_mass

print("x,y,z: %f,%f,%f" % (x,y,z))
cmd.pseudoatom("COM", pos=[x,y,z], name="COM", label="Center of mass")
cmd.show("spheres", "COM")
cmd.set("label_position", (0,3,0), "COM")
cmd.color("red", "COM")
cmd.set("sphere_scale", "1.01", "COM")
