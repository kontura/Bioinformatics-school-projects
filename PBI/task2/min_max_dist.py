from pymol import *

cmd.select("s", "all")
model = cmd.get_model("s", "1")
atoms = model.atom

maxi = 0.0
mini = dist=math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(atoms[0].coord,model.atom[2].coord))))

alls = model.atom
cutted = model.atom

for c1 in range(len(atoms)):
    a = atoms[c1]
    cmd.select("neighborhood", "neighbor id %s"%a.id)
    if (cmd.count_atoms("neighborhood") > 0):
        neighborhood_model = cmd.get_model("neighborhood", "1")
        neighborhood_atoms = neighborhood_model.atom
        neighborhood_atoms_ids = [elem.id for elem in neighborhood_atoms]
    else:
        neighborhood_atoms_ids = []

    for c2 in range(c1+1, len(atoms)):
        b = atoms[c2]
        dist=math.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(a.coord,b.coord))))
        if dist > maxi:
            maxi = dist
            maxi_one_id = a.id
            maxi_two_id = b.id
        if dist < mini and (b.id not in neighborhood_atoms_ids):
            mini = dist
            mini_one_id = a.id
            mini_two_id = b.id


cmd.select("one", "id %s"%maxi_one_id)
cmd.select("two", "id %s"%maxi_two_id)
cmd.color("red", "one")
cmd.color("red", "two")
cmd.distance("maximum", "one", "two")
print("maximum: %f" % maxi)

cmd.select("one", "id %s"%mini_one_id)
cmd.select("two", "id %s"%mini_two_id)
cmd.distance("minimum", "one", "two")
cmd.color("green", "one")
cmd.color("green", "two")
print("minimum: %f" % mini)
#cmd.show("spheres", "one")
#cmd.show("spheres", "two")
#cmd.set("sphere_scale", "0.01", "one")
#cmd.set("sphere_scale", "0.01", "two")
