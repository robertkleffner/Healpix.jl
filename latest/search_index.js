var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "index.html#Healpix.jl:-an-implementation-of-the-Healpix-tessellation-scheme-in-Julia-1",
    "page": "Introduction",
    "title": "Healpix.jl: an implementation of the Healpix tessellation scheme in Julia",
    "category": "section",
    "text": "This is the documentation of the Healpix.jl package, an implementation of the Healpix spherical tessellation scheme written entirely in Julia.This library is a work-in-progress: if you want something with more functionality, have a look at Libhealpix.jl, as it wraps the Healpix C++ library. This package has the main purpose of providing a Julia-only solution, so that it can easily be used on platforms not supported by the Healpix C++ library (e.g., Windows).This library implements algorithms for converting directions into pixel indices and vice versa. It supports both RING and NESTED schemes, and it employs Julia\'s powerful type system to avoid mistaking one scheme in place of the other."
},

{
    "location": "index.html#Documentation-1",
    "page": "Introduction",
    "title": "Documentation",
    "category": "section",
    "text": "The documentation was built using Documenter.jl.println(\"Documentation built $(now()) with Julia $(VERSION).\") # hide"
},

{
    "location": "index.html#Index-1",
    "page": "Introduction",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "resolutions.html#",
    "page": "Working with resolutions",
    "title": "Working with resolutions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "resolutions.html#Healpix.Resolution",
    "page": "Working with resolutions",
    "title": "Healpix.Resolution",
    "category": "type",
    "text": "Resolution objects are needed to perform a number of pixel-related functions, e.g., convert a direction into a pixel number and vice versa.\n\n\n\n"
},

{
    "location": "resolutions.html#Healpix.Resolution-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.Resolution",
    "category": "method",
    "text": "Resolution(nside) -> Resolution\n\nCreate a Resolution object, given a value for NSIDE.\n\n\n\n"
},

{
    "location": "resolutions.html#Healpix.nsideok-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nsideok",
    "category": "method",
    "text": "nsideok(nside::Integer) -> Bool\n\nCheck whether nside is a valid NSIDE parameter.\n\n\n\n"
},

{
    "location": "resolutions.html#Healpix.nside2npix-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nside2npix",
    "category": "method",
    "text": "nside2npix(nside::Integer) -> Integer\n\nReturn the number of pixels for a Healpix map with the specified NSIDE value. If NSIDE is not an integer power of two, the function throws a DomainError exception.\n\n\n\n"
},

{
    "location": "resolutions.html#Healpix.npix2nside-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.npix2nside",
    "category": "method",
    "text": "npix2nside(npix::Integer) -> Integer\n\nGiven the number of pixels in a Healpix map, return the NSIDE resolution parameter. If the number is invalid, throw a DomainError exception.\n\n\n\n"
},

{
    "location": "resolutions.html#Healpix.nside2pixarea-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nside2pixarea",
    "category": "method",
    "text": "nside2pixarea(nside::Integer) -> Real\n\nReturn the solid angle of a pixel in a map with the specified NSIDE parameter. The result is expressed in steradians.\n\n\n\n"
},

{
    "location": "resolutions.html#Healpix.nside2resol-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nside2resol",
    "category": "method",
    "text": "nside2resol(nside::Integer) -> Real\n\nReturn the approximate resolution of a map with the specified NSIDE. The resolution is expressed in radians, and it is the square root of the pixel size.\n\n\n\n"
},

{
    "location": "resolutions.html#Working-with-resolutions-1",
    "page": "Working with resolutions",
    "title": "Working with resolutions",
    "category": "section",
    "text": "A Healpix tessellation is parametrized by a number, called NSIDE, which must be a positive power of 2. It is related to the number of pixels N in the maps by the simple equation N = 12 mathrmNSIDE^2, and it is therefore related to the resolution of the pixelization. Any function working on a Healpix tessellation needs to receive the value of NSIDE. Healpix.jl provides a wrapper around this parameter, the Resolution type, which internally keeps a number of precomputed coefficients to accelerate calculations.The following example prints a table containing details about a few Healpix resolutions:using Healpix # hide\n@printf(\"%-6s\\t%-12s\\t%-12s\\t%-12s\\n\",\n        \"NSIDE\",\n        \"#pix\",\n        \"#pix per face\",\n        \"solid angle\")\nfor poweroftwo in [0, 1, 2, 3, 4, 5]\n    res = Resolution(2 ^ poweroftwo)\n    @printf(\"%6d\\t%12d\\t%12d\\t%12.4f\\n\",\n            res.nside,\n            res.numOfPixels,\n            res.pixelsPerFace,\n            4π / res.numOfPixels)\nendResolution\nResolution(nside::Integer)\nnsideok(nside::Integer)\nnside2npix(nside::Integer)\nnpix2nside(npix::Integer)\nnside2pixarea(nside::Integer)\nnside2resol(nside::Integer)"
},

{
    "location": "pixelfunc.html#",
    "page": "Pixel functions",
    "title": "Pixel functions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "pixelfunc.html#Healpix.ang2vec-Tuple{Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.ang2vec",
    "category": "method",
    "text": "ang2vec(theta, phi) -> (Float64, Float64, Float64)\n\nGiven a direction in the sky with colatitude theta and longitude phi (in radians), return a tuple containing the x, y, and z components of the one-length vector pointing to that direction.\n\n\n\n"
},

{
    "location": "pixelfunc.html#Healpix.vec2ang-Tuple{Any,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.vec2ang",
    "category": "method",
    "text": "vec2ang(x, y, z) -> (Float64, Float64)\n\nGiven a vector (not necessarily normalized) whose Cartesian components are x, y, and z, return a pair (theta, phi) containing the colatitude theta and the longitude phi (in radians) of the direction in the sky the vector is pointing at.\n\n\n\n"
},

{
    "location": "pixelfunc.html#Healpix.ang2pixNest-Tuple{Healpix.Resolution,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.ang2pixNest",
    "category": "method",
    "text": "ang2pixNest(resol::Resolution, theta, phi) -> Integer\n\nReturn the index of the pixel which contains the point with coordinates (theta, the colatitude, and phi, the longitude), in radians, for a Healpix map with pixels in nested order. Note that pixel indexes are 1-based (this is Julia)!\n\n\n\n"
},

{
    "location": "pixelfunc.html#Healpix.ang2pixRing-Tuple{Healpix.Resolution,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.ang2pixRing",
    "category": "method",
    "text": "ang2pixRing(resol::Resolution, theta, phi) -> Integer\n\nReturn the index of the pixel which contains the point with coordinates (theta, the colatitude, and phi, the longitude), in radians, for a Healpix map with pixels in ring order. Note that pixel indexes are 1-based (this is Julia)!\n\n\n\n"
},

{
    "location": "pixelfunc.html#Healpix.pix2angNest-Tuple{Healpix.Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.pix2angNest",
    "category": "method",
    "text": "pix2angNest(resol::Resolution, pixel) -> (Float64, Float64)\n\nGiven the (1-based) index of a pixel in a Healpix map in nested order, return a pair containing the (colatitude, longitude) angles corresponding to its center, both expressed in radians.\n\n\n\n"
},

{
    "location": "pixelfunc.html#Healpix.pix2angRing-Tuple{Healpix.Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.pix2angRing",
    "category": "method",
    "text": "pix2angRing(resol::Resolution, pixel) -> (Float64, Float64)\n\nGiven the (1-based) index of a pixel in a Healpix map in ring order, return a pair containing the (colatitude, longitude) angles corresponding to its center, both expressed in radians.\n\n\n\n"
},

{
    "location": "pixelfunc.html#Pixel-functions-1",
    "page": "Pixel functions",
    "title": "Pixel functions",
    "category": "section",
    "text": "In this section we document the functions that convert from a direction in the sky into a pixel index, and vice versa.First of all, Healpix.jl implements the most basic functions to convert between spherical and Cartesian coordinates. Note that Healpix uses co-latitude instead of latitude:using Healpix # hide\nang2vec(0.0, 0.0)\nvec2ang(0.0, 0.0, 1.0)More interesting functions return the index of the pixel on a Healpix-tessellated sphere. For these functions to work, you have to provide a Resolution object:res = Resolution(16)\nang2pixRing(res, π/2, 0)\nang2pixNest(res, π/2, 0)ang2vec(theta, phi)\nvec2ang(x, y, z)\nang2pixNest(resol::Resolution, theta, phi)\nang2pixRing(resol::Resolution, theta, phi)\npix2angNest(resol::Resolution, pixel)\npix2angRing(resol::Resolution, pixel)"
},

{
    "location": "mapfunc.html#",
    "page": "Map functions",
    "title": "Map functions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "mapfunc.html#Healpix.Map",
    "page": "Map functions",
    "title": "Healpix.Map",
    "category": "type",
    "text": "A Healpix map. The type T is used for the value of the pixels in a map, and it can be anything (even a string!). The type O is used to specify the ordering of the pixels, and it can either be RingOrder or NestedOrder.\n\n\n\n"
},

{
    "location": "mapfunc.html#Map-functions-1",
    "page": "Map functions",
    "title": "Map functions",
    "category": "section",
    "text": "Functions like pix2angNest and ang2pixNest fully define the Healpix tessellation scheme. They are however extremely impractical in a number of situations. It happens often that a large fraction of pixels in a map need to be processed together. Healpix.jl introduces the Map{T, O <: Order} type, which acts as a collection of all the pixels on the sphere. A Map type holds the value of all the pixels in its pixels field, and it keeps track of the ordering (either RING or NESTED). Here is an example that shows how to create a map and initialize it:nside = 32\nm = Map{Float64, RingOrder}(nside)\nm.pixels[:] = 1.0  # Set all pixels to 1Map"
},

{
    "location": "mapfunc.html#Healpix.Order",
    "page": "Map functions",
    "title": "Healpix.Order",
    "category": "type",
    "text": "Abstract type representing the ordering of pixels in a Healpix map. See also RingOrder and NestedOrder.\n\n\n\n"
},

{
    "location": "mapfunc.html#Healpix.RingOrder",
    "page": "Map functions",
    "title": "Healpix.RingOrder",
    "category": "type",
    "text": "The RingOrder type should be used when creating Map types in order to specify that the pixels in the map are sorted in ring ordering. (See also NestedOrder.)\n\n\n\n"
},

{
    "location": "mapfunc.html#Healpix.NestedOrder",
    "page": "Map functions",
    "title": "Healpix.NestedOrder",
    "category": "type",
    "text": "The NestedOrder type should be used when creating Map types in order to specify that the pixels in the map are sorted in ring ordering. (See also RingOrder.)\n\n\n\n"
},

{
    "location": "mapfunc.html#Encoding-the-order-1",
    "page": "Map functions",
    "title": "Encoding the order",
    "category": "section",
    "text": "Healpix.jl distinguishes between RING and NEST orderings using Julia\'s typesystem. The abstract type Order has two descendeants, RingOrder and NestedOrder, which are used to instantiate objects of type Map.Order\nRingOrder\nNestedOrder"
},

{
    "location": "mapfunc.html#Healpix.pix2ang",
    "page": "Map functions",
    "title": "Healpix.pix2ang",
    "category": "function",
    "text": "pix2ang{T, O <: Order}(map::Map{T, O}, ipix) -> (Float64, Float64)\n\nReturn the pair (theta, phi), where theta is the colatitude and phi the longitude of the direction of the pixel center with index ipix.\n\n\n\n"
},

{
    "location": "mapfunc.html#Healpix.ang2pix",
    "page": "Map functions",
    "title": "Healpix.ang2pix",
    "category": "function",
    "text": "ang2pix{T, O <: Order}(map::Map{T, O}, theta::Real, phi::Real)\n\nConvert the direction specified by the colatitude theta (∈ [0, π]) and the longitude phi (∈ [0, 2π]) into the index of the pixel in the Healpix map map.\n\n\n\n"
},

{
    "location": "mapfunc.html#Pixel-functions-1",
    "page": "Map functions",
    "title": "Pixel functions",
    "category": "section",
    "text": "When working with maps, it is not needed to pick between ang2pixNest and ang2pixRing because a Map type already encodes the ordering. Functions pix2ang and ang2pix always choose the correct ordering, but they require a Map instead of a Resolution as their first argument.pix2ang\nang2pix"
},

{
    "location": "mapfunc.html#Healpix.saveToFITS",
    "page": "Map functions",
    "title": "Healpix.saveToFITS",
    "category": "function",
    "text": "saveToFITS{T <: Number, O <: Order}(map::Map{T, O},\n                                    f::FITSIO.FITSFile,\n                                    column)\nsaveToFITS{T <: Number, O <: Order}(map::Map{T, O},\n                                    fileName::String,\n                                    typechar=\"D\",\n                                    unit=\"\",\n                                    extname=\"MAP\")\n\nSave a Healpix map in the specified (1-based index) column in a FITS file. If the code fails, FITSIO will raise an exception. (Refer to the FITSIO library for more information.)\n\n\n\nsaveToFITS(map::Map{T, O}, filename::AbstractString, typechar=\"D\", unit=\"\", extname=\"MAP\") where {T <: Number, O <: Order}\n\nSave a map into a FITS file. The name of the file is specified in filename; if it begins with !, existing files will be overwritten without warning. The parameter typechar specifies the data type to be used in the FITS file: the default (D) will save 64-bit floating-point values. See the CFITSIO documentation for other values. The keyword unit specifies the measure unit used for the pixels in the map. The keyword extname specifies the name of the HDU where the map pixels will be written.\n\n\n\n"
},

{
    "location": "mapfunc.html#Healpix.savePixelsToFITS",
    "page": "Map functions",
    "title": "Healpix.savePixelsToFITS",
    "category": "function",
    "text": "savePixelsToFITS(map::Map{T}, f::FITSIO.FITSFile, column) where {T <: Number}\n\nSave the pixels of map into the column with index/name column in the FITS file, which must have been already opened.\n\n\n\n"
},

{
    "location": "mapfunc.html#Healpix.readMapFromFITS",
    "page": "Map functions",
    "title": "Healpix.readMapFromFITS",
    "category": "function",
    "text": "readMapFromFITS{T <: Number}(f::FITSIO.FITSFILE, column, t::Type{T})\nreadMapFromFITS{T <: Number}(fileName::String, column, t::Type{T})\n\nRead a Healpix map from the specified (1-base indexed) column in a FITS file. The values will be read as numbers of type T. If the code fails, FITSIO will raise an exception. (Refer to the FITSIO library for more information.)\n\n\n\n"
},

{
    "location": "mapfunc.html#Loading-and-saving-maps-1",
    "page": "Map functions",
    "title": "Loading and saving maps",
    "category": "section",
    "text": "Healpix.jl implements a number of functions to save maps in FITS files.saveToFITSFunction savePixelsToFITS is a low-level function. It knows nothing about the ordering schema used for the pixels, so the caller should manually write the ORDERING keyword in the HDU header by itself.savePixelsToFITSTo load a map from a FITS file, you can use readMapFromFITS.readMapFromFITS"
},

{
    "location": "mapfunc.html#Healpix.conformables",
    "page": "Map functions",
    "title": "Healpix.conformables",
    "category": "function",
    "text": "conformables{T, S, O1 <: Order, O2 <: Order}(map1::Map{T, O1},\n                                             map2::Map{S, O2}) -> Bool\n\nDetermine if two Healpix maps are \"conformables\", i.e., if their shape and ordering are the same.\n\n\n\n"
},

{
    "location": "mapfunc.html#Testing-for-conformability-1",
    "page": "Map functions",
    "title": "Testing for conformability",
    "category": "section",
    "text": "It often happens that two Healpix maps need to be combined together: for instance, pixels on a sky map might need to be masked using a sky mask, or one map might need to be subtracted from another one. «Conformability» means that the operation between the two maps can be done directly on the pixels, without oordering or resolution conversions. The function conformables checks this.m1 = Map{Float64, RingOrder}(1)\nm2 = Map{Float64, RingOrder}(1)\nm3 = Map{Float64, NestedOrder}(1)\nm4 = Map{Float64, NestedOrder}(2)\nconformables(m1, m2)\nconformables(m1, m3)\nconformables(m1, m4)conformables"
},

{
    "location": "alm.html#",
    "page": "Spherical harmonics",
    "title": "Spherical harmonics",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "alm.html#Healpix.Alm",
    "page": "Spherical harmonics",
    "title": "Healpix.Alm",
    "category": "type",
    "text": "An array of a_ℓm numbers.\n\n\n\n"
},

{
    "location": "alm.html#Healpix.numberOfAlms",
    "page": "Spherical harmonics",
    "title": "Healpix.numberOfAlms",
    "category": "function",
    "text": "numberOfAlms(lmax::Integer, mmax::Integer) -> Integer\nnumberOfAlms(lmax::Integer) -> Integer\n\nReturn the size of the array of complex numbers needed to store the a_lm coefficients in the range of ℓ and m specified by lmax and mmax. If mmax is not specified, it is assumed to be equal to lmax. If lmax and mmax are inconsistent or negative, a DomainError exception is thrown.\n\n\n\n"
},

{
    "location": "alm.html#Spherical-harmonics-1",
    "page": "Spherical harmonics",
    "title": "Spherical harmonics",
    "category": "section",
    "text": "The support for spherical harmonics in Healpix.jl is still woefully inadequate. Only a few functions to load and store harmonic coefficients are available. Everything revolves around the Alm type:AlmThe number of coefficients in a spherical harmonic expansion is infinite. For obvious reasons, Healpix.jl only allows to store band-limited expansions. The function numberOfAlms returns the number of floating-point numbers used to store the expansion, as a function of the maximum value for ell and m.numberOfAlms"
},

{
    "location": "alm.html#Healpix.readAlmFromFITS",
    "page": "Spherical harmonics",
    "title": "Healpix.readAlmFromFITS",
    "category": "function",
    "text": "readAlmFromFITS{T <: Complex}(f::FITSIO.FITSFile, t::Type{T}) -> Alm{T}\nreadAlmFromFITS{T <: Complex}(fileName::String, t::Type{T}) -> Alm{T}\n\nRead a set of a_ℓm coefficients from a FITS file. If the code fails, FITSIO will raise an exception. (Refer to the FITSIO library for more information.)\n\n\n\n"
},

{
    "location": "alm.html#Loading-and-saving-harmonic-coefficients-1",
    "page": "Spherical harmonics",
    "title": "Loading and saving harmonic coefficients",
    "category": "section",
    "text": "readAlmFromFITS"
},

{
    "location": "visualization.html#",
    "page": "Visualization",
    "title": "Visualization",
    "category": "page",
    "text": ""
},

{
    "location": "visualization.html#Healpix.project",
    "page": "Visualization",
    "title": "Healpix.project",
    "category": "function",
    "text": "project(m::Map{T, O}; kwargs...) where {T, O <: Order}\n\nReturn a 2D bitmap (array) containing a cartographic projection of the map and a 2D bitmap containing a boolean mask. The size of the bitmap is specified by figsize, which must be a 2-tuple. The function projfn must be a function which accepts as input two parameters x and y (numbers between -1 and 1).\n\nThe following keywords can be used in the call:\n\nfigsize: 2-tuple specifying the (height, width) of the bitmap in pixels\ncenter: 2-tuple specifying the location (colatitude, longitude) of the sky point that is to be placed in the middle of the image (in radians)\nshow: Boolean; if true (the default), the bitmap will be shown using functions from the \"Plots\" package.\nreturnmask: Boolean; if true, the function returns a 2-tuple containing the image bitmap and a mask bitmap, which is set to true if the pixel falls within the carthographic projection, false otherwise. If returnmask is false, only the image bitmap is returned.\nshow: Boolean. If true (the default), the map will be displayed. It has no effect if returnmask is true.  \n\n\n\n"
},

{
    "location": "visualization.html#Visualization-functions-1",
    "page": "Visualization",
    "title": "Visualization functions",
    "category": "section",
    "text": "Healpix.jl implements the project function, which creates a 2D matrix containing the cartographic projection of a map. A few standard cartographic projections are implemented, but users can provide their own projections. The function can optionally use heatmap (from the Plots.jl package) to display the 2D matrix. Two useful wrappers to project are equirectangular and mollweide, which employ the equirectangular and Mollweide projections respectively.using Healpix\n\nnside = 8\nm = Map{Float64, RingOrder}(nside)\nm.pixels[:] = 1:length(m.pixels)\nmollweide(m)project"
},

{
    "location": "visualization.html#Healpix.equiprojinv",
    "page": "Visualization",
    "title": "Healpix.equiprojinv",
    "category": "function",
    "text": "function equiprojinv(x, y)\n\nInverse equirectangular projection. Given a point (x, y) on the plane [-1, 1] × [-1, 1], return a tuple (Bool, Number, Number) where the first Boolean is a flag telling if the point falls within the projection (true) or not (false), and the two numbers are the latitude and colatitude in radians.\n\n\n\n"
},

{
    "location": "visualization.html#Healpix.mollweideprojinv",
    "page": "Visualization",
    "title": "Healpix.mollweideprojinv",
    "category": "function",
    "text": "function mollweideprojinv(x, y)\n\nInverse Mollweide projection. Given a point (x, y) on the plane, with x ∈ [-1, 1], y ∈ [-1, 1], return a 3-tuple of type (Bool, Number, Number). The boolean specifies if (x, y) falls within the map (true) or not (false), the second and third arguments are the latitude and longitude in radians.\n\n\n\n"
},

{
    "location": "visualization.html#Cartographic-projections-1",
    "page": "Visualization",
    "title": "Cartographic projections",
    "category": "section",
    "text": "equiprojinv\nmollweideprojinv"
},

{
    "location": "visualization.html#Healpix.equirectangular",
    "page": "Visualization",
    "title": "Healpix.equirectangular",
    "category": "function",
    "text": "equirectangular(m::Map{T,O}; kwargs...) where {T <: AbstractFloat, O <: Order}\n\nHigh-level wrapper around project for equirectangular projections.\n\n\n\n"
},

{
    "location": "visualization.html#Healpix.mollweide",
    "page": "Visualization",
    "title": "Healpix.mollweide",
    "category": "function",
    "text": "mollweide(m::Map{T,O}; kwargs...) where {T <: AbstractFloat, O <: Order}\n\nHigh-level wrapper around project for Mollweide projections.\n\n\n\n"
},

{
    "location": "visualization.html#High-level-wrappers-to-project-1",
    "page": "Visualization",
    "title": "High-level wrappers to project",
    "category": "section",
    "text": "equirectangular\nmollweide"
},

{
    "location": "misc.html#",
    "page": "Miscellanea",
    "title": "Miscellanea",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "misc.html#Healpix.normalizeAngle-Tuple{Any}",
    "page": "Miscellanea",
    "title": "Healpix.normalizeAngle",
    "category": "method",
    "text": "normalizeAngle(x)\n\nReturn the same angle as the argument, but in the range [0, 2π). Note that this is slightly different from mod2pi, as the latter returns a value in the range [0, 2π].\n\n\n\n"
},

{
    "location": "misc.html#Healpix.lat2colat-Tuple{Any}",
    "page": "Miscellanea",
    "title": "Healpix.lat2colat",
    "category": "method",
    "text": "lat2colat(x)\n\nConvert latitude into colatitude. Both x and the result are expressed in radians.\n\n\n\n"
},

{
    "location": "misc.html#General-purpose-functions-1",
    "page": "Miscellanea",
    "title": "General-purpose functions",
    "category": "section",
    "text": "Healpix.jl implements a few generic functions that can be helpful when doing calculations on the sphere.normalizeAngle(x)\nlat2colat(x)"
},

]}
