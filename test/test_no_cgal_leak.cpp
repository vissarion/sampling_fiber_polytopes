// Structural guard: no public header may leak CGAL / LEDA / respol.
//
// The PIMPL boundary in include/secpoly/opt_oracle.h is what allows the rest of
// the project -- and volesti -- to compile without ever seeing respol's patched
// CGAL Kernel_d. If a public header ever pulls CGAL in, two different definitions
// of the same CGAL classes (patched vs stock) could end up in one program, which
// is silent UB rather than a build error. So we fail loudly here instead.
//
// NOTE: this TU is deliberately compiled WITHOUT the respol include directories.

#include <secpoly/opt_oracle.h>

#if defined(CGAL_VERSION) || defined(CGAL_VERSION_NR) || defined(CGAL_BASIC_H)
#error "CGAL leaked into a public secpoly header -- the PIMPL boundary is broken"
#endif

#if defined(LEDA_ROOT) || defined(LEDA_BEGIN_NAMESPACE)
#error "LEDA leaked into a public secpoly header -- the PIMPL boundary is broken"
#endif

#include <cstdio>

int main() {
    std::printf("test_no_cgal_leak: PASS\n");
    return 0;
}
