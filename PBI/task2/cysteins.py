from pymol import *

cmd.color("tv_blue", "all")

cmd.select("cysteins", "resn Cys")
cmd.color("gold", "cysteins")

model = cmd.get_model("cysteins", "1")

for source in model.atom:
    cmd.select("one", "id %s"%source.id)
    for target in model.atom:
        cmd.select("two", "id %s"%target.id)
        dist = cmd.get_distance("one","two")
        if (dist > 1.01 and dist < 3):
            cmd.distance("d", "one", "two")
            cmd.color("orange", "one")
            cmd.color("orange", "two")
            cmd.label("one", "index")
            cmd.label("two", "index")
            cmd.hide("dashes", "d")
            cmd.set("label_color", "orange", "d")
