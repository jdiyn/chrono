// SWIG forwards byte literals into uint expressions and the wrapper build fails (CS0266)
// Override with explicit byte casts fixes this:
%feature("cs:constvalue", "(byte)0") chrono::UpdateFlags::NONE;
%feature("cs:constvalue", "(byte)(1u << 0)") chrono::UpdateFlags::DYNAMICS;
%feature("cs:constvalue", "(byte)(1u << 1)") chrono::UpdateFlags::JACOBIANS;
%feature("cs:constvalue", "(byte)(1u << 2)") chrono::UpdateFlags::VISUAL_ASSETS;
%feature("cs:constvalue", "(byte)((1u << 0) | (1u << 1) | (1u << 2))") chrono::UpdateFlags::UPDATE_ALL;
%feature("cs:constvalue", "(byte)((1u << 0) | (1u << 1))") chrono::UpdateFlags::UPDATE_ALL_NO_VISUAL;

// UpdateFlags must be known before parsing headers (e.g., ChMesh.h) that use it.
%include "../../../chrono/physics/ChUpdateFlags.h"

// explicit alias for fea so UpdateFlags resolves in that namespace (sometimes threw errors otherwise)
%inline %{
namespace chrono {
namespace fea {
typedef chrono::UpdateFlags UpdateFlags;
}
}
%}
