import kroot.toms748
import kotlin.math.*
import kotlin.test.Test
import kotlin.test.assertTrue

class TOMS748 {
    @Test
    fun f1() {
        val root = toms748(PI/2, PI, 1E-15, 100) { x -> sin(x) - x / 2 }
        assertTrue("Incorrect value") {
            abs(root.x - 1.89549) < 1E-5
        }
        assertTrue("Too many iterations") {
            root.iterations <= 10
        }
    }

    @Test
    fun f3() {
        val a = -100
        val b = -2
        val root = toms748(-9.0, 31.0, 1E-15, 100) { x -> a * x * exp(b * x) }
        assertTrue("Incorrect value") {
            abs(root.x) < 1E-5
        }
        assertTrue("Too many iterations") {
            root.iterations <= 26
        }
    }

    @Test
    fun f12() {
        val root = toms748(1.0, 100.0, 1E-15, 100) { x -> sqrt(x) - sqrt(2.0) }
        assertTrue("Incorrect value") {
            abs(root.x - 2) < 1E-5
        }
        assertTrue("Too many iterations") {
            root.iterations <= 5
        }
    }
}