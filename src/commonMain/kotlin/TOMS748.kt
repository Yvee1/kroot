@file:Suppress("NOTHING_TO_INLINE")
package kroot

import kotlin.math.abs
import kotlin.math.sign

const val EPS = 1E-14

data class Root(val x: Double, val iterations: Int)

/**
 * Returns the approximate root of [f] within interval [[a], [b]].
 *
 * Assumes [f] is continuous and ```f(a) * f(b) < 0```.
 * Implements Algorithm 748 by Alefeld, Potra and Shi in Transactions on Mathematical Software (TOMS).
 */
fun toms748(a: Double, b: Double, tol: Double, maxIter: Int, f: (Double) -> Double) =
    toms748(a, b, f(a), f(b), tol, maxIter, f)

fun toms748(a: Double, b: Double, fa: Double, fb: Double, tol: Double, maxIter: Int, f: (Double) -> Double): Root {
    require(a < b) {
        "Parameter a should be strictly smaller than b.\na: $a\nb: $b"
    }
    require(tol > 0) {
        "The tolerance should be positive.\ntol: $tol"
    }
    require(maxIter > 0) {
        "The maximum number of iterations should be positive.\nmaxIter: $maxIter"
    }

    var iters = 0

    fun root(x: Double) = Root(x, iters)

    if (fa == 0.0) return root(a)
    if (fb == 0.0) return root(b)
    if (a >= b - tol) return root((a + b) / 2)

    require(sign(fa) * sign(fb) < 0) {
        "Parameters a and b do not bracket the root: their function values should have opposite signs.\nfa: $fa\nfb: $fb"
    }

    // First step uses the secant method
    var ax = a
    var bx = b
    var fax = fa
    var fbx = fb
    var cx = secant(a, b, fax, fbx)
    var fcx = f(cx)
    iters++

    // Return if done
    if (fcx == 0.0 || iters >= maxIter) return root(cx)

    // Update bracket
    val bracket = bracket(ax, bx, cx, fax, fbx, fcx)
    ax = bracket.a
    bx = bracket.b
    var dx = bracket.d
    fax = bracket.fa
    fbx = bracket.fb
    var fdx = bracket.fd

    // Second step uses quadratic interpolation
    cx = newtonQuadratic(ax, bx, dx, fax, fbx, fdx, 2)
    iters++
    var ex = dx
    var fex = fdx

    // Return if done
    fcx = f(cx)
    if (fcx == 0.0 || iters >= maxIter) return root(cx)

    // Update bracket
    val brack = bracket(ax, bx, cx, fax, fbx, fcx)
    ax = brack.a
    bx = brack.b
    dx = brack.d
    fax = brack.fa
    fbx = brack.fb
    fdx = brack.fd

    while (iters < maxIter && ax <= bx - tol) {
        val a0 = ax
        val b0 = bx

        for (nSteps in 2..3) {
            // Cubic interpolation if the four function values are distinct. Otherwise, a quadratic step is taken.
            cx = if (distinct(fax, fbx, fdx, fex)) {
                val r = ipzero(ax, bx, dx, ex, fax, fbx, fdx, fex)
                if (ax < r && r < bx) {
                    r
                } else {
                    newtonQuadratic(ax, bx, dx, fax, fbx, fdx, nSteps)
                }
            } else {
                newtonQuadratic(ax, bx, dx, fax, fbx, fdx, nSteps)
            }

            iters++
            fcx = f(cx)

            // Return if done
            if (fcx == 0.0 || iters >= maxIter) return root(cx)

            // Update bracket
            ex = dx
            fex = fdx

            val br = bracket(ax, bx, cx, fax, fbx, fcx)
            ax = br.a
            bx = br.b
            dx = br.d
            fax = br.fa
            fbx = br.fb
            fdx = br.fd
        }

        // Double-length secant step
        val (ux, fux) = if (abs(fax) < abs(fbx)) ax to fax else bx to fbx
        cx = ux - 2 * fux / f2(ax, bx, fax, fbx)
        if (abs(cx - ux) > (bx - ax)/2) {
            cx = (bx + ax) / 2
        }

        iters++
        fcx = f(cx)

        // Return if done
        if (fcx == 0.0 || iters >= maxIter) return root(cx)

        // Update bracket
        val br = bracket(ax, bx, cx, fax, fbx, fcx)
        ax = br.a
        bx = br.b
        dx = br.d
        fax = br.fa
        fbx = br.fb
        fdx = br.fd

        // Are we converging?
        if (bx - ax < (b0 - a0) / 2)
            continue

        // If not, take a bisection step
        ex = dx
        fex = fdx
        val z = (a + b) / 2
        iters++
        val fz = f(z)
        if (fz == 0.0 || iters >= maxIter) return root(z)
        val brr = bracket(ax, bx, z, fax, fbx, fz)
        ax = brr.a
        bx = brr.b
        dx = brr.d
        fax = brr.fa
        fbx = brr.fb
        fdx = brr.fd
    }

    return root((ax + bx)/2)
}

private inline fun distinct(fa: Double, fb: Double, fd: Double, fe: Double): Boolean =
    abs(fa - fb) > EPS && abs(fa - fd) > EPS && abs(fa - fe) > EPS &&
            abs(fb - fd) > EPS && abs(fb - fe) > EPS && abs(fd - fe) > EPS

data class Bracket(val a: Double, val b: Double, val d: Double, val fa: Double, val fb: Double, val fd: Double)

private inline fun bracket(a: Double, b: Double, c: Double, fa: Double, fb: Double, fc: Double): Bracket =
    if (sign(fa) * sign(fc) < 0) Bracket(a, c, b, fa, fc, fb)
    else Bracket(c, b, a, fc, fb, fa)

/**
 * Returns ```f[a, b] := (f(b) - f(a))/(b-a)```, the mean of [f] on interval [[a], [b]]
 */
private inline fun f2(a: Double, b: Double, fa: Double, fb: Double): Double =
    (fb - fa) / (b - a)

/**
 * Returns ```f[a, b, d] := (f[b, d] - f[a, b])/(d-a)```.
 */
private inline fun f3(a: Double, b: Double, d: Double, fa: Double, fb: Double, fd: Double) =
    (f2(b, d, fb, fd) - f2(a, b, fa, fb)) / (d - a)

/**
 * Returns the unique zero of the quadratic polynomial
 * ```f[a, b, d](x-a)(x-b) + f[a, b](x-a) + f(a)```.
 */
private fun newtonQuadratic(a: Double, b: Double, d: Double, fa: Double, fb: Double, fd: Double, k: Int): Double {
    val A = f3(a, b, d, fa, fb, fd)
    val B = f2(a, b, fa, fb)

    return if (A == 0.0) a - fa / B
    else {
        var r = if (A * fa > 0) a else b
        for (i in 1..k) {
            val rn = r - ((A * (r - b) + B) * (r - a) + fa) / (B + A * (2 * (r - a - b)))
            r = if (a < rn && rn < b) {
                rn
            } else {
                if (a < r && r < b)
                    return r
                (a + b) / 2
            }
        }
        r
    }
}

private fun ipzero(a: Double, b: Double, c: Double, d: Double, fa: Double, fb: Double, fc: Double, fd: Double): Double {
    val q11 = (c - d) * fc / (fd - fc)
    val q21 = (b - c) * fb / (fc - fb)
    val q31 = (a - b) * fa / (fb - fa)
    val d21 = (b - c) * fc / (fc - fb)
    val d31 = (a - b) * fb / (fb - fa)
    val q22 = (d21 - q11) * fb / (fd - fb)
    val q32 = (d31 - q21) * fa / (fc - fa)
    val d32 = (d31 - q21) * fc / (fc - fa)
    val q33 = (d32 - q22) * fa / (fd - fa)

    return a + q31 + q32 + q33
}

private inline fun secant(a: Double, b: Double, fa: Double, fb: Double): Double {
    val c = a - fa / f2(a, b, fa, fb)
    return if (c <= a + abs(a) * EPS || c >= b - abs(b) * EPS) (a + b) / 2 else c
}