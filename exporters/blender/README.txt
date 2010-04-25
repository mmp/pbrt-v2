This is a very alpha-quality exporter from blender for the latest (alpha 2.0)
version of the pbrt renderer (which is available from http://github.com/mmp/pbrt-v2).

Current TODOs/bugs:
- testing!  Many materials, integrators, etc, haven't been well-tested
- bug: output image filename is "default", not "default.exr"
- pbrt v2 has some built-in data with subsurface scattering properties for various
  media.  need to get these wired up to the GUI
- similarly the way that pbrt gets IOR and k values for the metal material is
  from text files; need to wire those up in the exporter

lux features that pbrt doesn't support:
- exr image output always
- no compositiong
- no portals
- no light groups
- autofocus, focus on object
- bokeh stuff
- camera lens shift
- image display interval, write interval
- color space
- tone mapping
- gamma correction
- no sunsky
- no ies lights
- no flipz
- only exr (and tga) texture files
- cloud volume

no interactive view of progress etc

(probably many more)

