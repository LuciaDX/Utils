<?php

namespace onyx\utils;

class SimplexNoise {
    public $grad3;
    public $p;
    public $perm;
    public function __construct(public $seed){
        if (!$seed) $seed = Mathf::Random();
        $this->grad3 = [[1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0],
            [1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1],
            [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1]];
        $this->p = [];
        for ($i = 0; $i < 256; $i++) {
            $this->p[$i] = floor($seed * 256);
        }
        $this->perm = [];
        for ($i = 0; $i < 512; $i++) {
            $this->perm[$i] = $this->p[$i & 255];
        }
    }

    public function Dot($g, $x, $y) {
        return $g[0]*$x + $g[1]*$y;
    }

    public function Noise($xin, $yin) {
        $F2 = 0.5 * (sqrt(3.0) - 1.0);
        $s = ($xin + $yin) * $F2;
        $i = floor($xin + $s);
        $j = floor($yin + $s);
        $G2 = (3.0 - sqrt(3.0)) / 6.0;
        $t = ($i + $j) * $G2;
        $X0 = $i - $t;
        $Y0 = $j - $t;
        $x0 = $xin - $X0;
        $y0 = $yin - $Y0;
        if ($x0 > $y0) {
            $i1 = 1;
            $j1 = 0;
        } else {
            $i1 = 0;
            $j1 = 1;
        }
        $x1 = $x0 - $i1 + $G2;
        $y1 = $y0 - $j1 + $G2;
        $x2 = $x0 - 1.0 + 2.0 * $G2;
        $y2 = $y0 - 1.0 + 2.0 * $G2;
        $ii = $i & 255;
        $jj = $j & 255;
        $gi0 = $this->perm[$ii + $this->perm[$jj]] % 12;
        $gi1 = $this->perm[$ii + $i1 + $this->perm[$jj + $j1]] % 12;
        $gi2 = $this->perm[$ii + 1 + $this->perm[$jj + 1]] % 12;
        $t0 = 0.5 - $x0 * $x0 - $y0 * $y0;
        if ($t0 < 0) {
            $n0 = 0.0;
        } else {
            $t0 *= $t0;
            $n0 = $t0 * $t0 * $this->Dot($this->grad3[$gi0], $x0, $y0);
        }
        $t1 = 0.5 - $x1 * $x1 - $y1 * $y1;
        if ($t1 < 0) {
            $n1 = 0.0;
        } else {
            $t1 *= $t1;
            $n1 = $t1 * $t1 * $this->Dot($this->grad3[$gi1], $x1, $y1);
        }
        $t2 = 0.5 - $x2 * $x2 - $y2 * $y2;
        if ($t2 < 0) {
            $n2 = 0.0;
        } else {
            $t2 *= $t2;
            $n2 = $t2 * $t2 * $this->Dot($this->grad3[$gi2], $x2, $y2);
        }
        return 70.0 * ($n0 + $n1 + $n2);
    }
}

class Mathf {

    public const FULL_ANGLE = 360;
    public const STRAIGHT_ANGLE = 180;
    public const GAMMA_TO_LINEAR = 2.2;
    public const LINEAR_TO_GAMMA = 0.45454545;
    public const IS_INTEGER = 0.5;
    public const RANDOM_SEED = 0.8694896071683615;

    public static function Epsilon(): float|int{
        return pow(2, -52);
    }

    public static function Random(): float|int{
        return (float)rand()/(float)getrandmax();
    }

    public static function Clamp($value, $min, $max): float|int{
        return $value < $min ? $min : ($value > $max ? $max : $value);
    }

    public static function InverseLerp($a, $b, $value): float|int{
        return (self::Clamp($value, min($a, $b), max($a, $b)) - $a) / ($b - $a);
    }

    public static function Approximately($f1, $f2): float|int{
        return abs($f1 - $f2) < self::Epsilon();
    }

    public static function ClosestPowerOfTwo($value): float|int{
        $nextPowerOfTwo = self::NextPowerOfTwo($value);
        if ($nextPowerOfTwo - $value > $nextPowerOfTwo >> 2) {
            return $nextPowerOfTwo >> 1;
        }
        return $nextPowerOfTwo;
    }

    public static function DeltaAngle($current, $target): float|int{
        if (abs($current) > self::FULL_ANGLE) {
            $current %= self::FULL_ANGLE;
        }
        if (abs($target) > self::FULL_ANGLE) {
            $target %= self::FULL_ANGLE;
        }
        return $target - $current;
    }

    public static function GammaToLinearSpace($value) {
        return pow($value, self::GAMMA_TO_LINEAR);
    }

    public static function Lerp($a, $b, $t) {
        return ($b - $a) * self::Clamp($t,0,1) + $a;
    }

    public static function LerpAngle($a, $b, $t) {
        while ($a > $b + self::STRAIGHT_ANGLE) {
            $b += self::FULL_ANGLE;
        }
        while ($b > $a + self::STRAIGHT_ANGLE) {
            $b -= self::FULL_ANGLE;
        }
        return self::Lerp($a, $b, $t);
    }

    public static function IsPowerOfTwo($value): float|int{
        return ($value & ($value - 1)) === 0;
    }

    public static function NextPowerOfTwo($value): float|int{
        if ($value < 0) return 0;
        --$value;
        $value |= $value >> 1;
        $value |= $value >> 2;
        $value |= $value >> 4;
        $value |= $value >> 8;
        $value |= $value >> 16;
        $value += 1;
        return $value;
    }

    public static function LerpUnclamped($a, $b, $t) {
        if ($t < 0 || $t > 1) {
            return $a + abs($b - $a) * $t;
        }
        return ($b - $a) * self::Clamp($t,0,1) + $a;
    }

    public static function LinearToGammaSpace($value) {
        return pow($value, self::LINEAR_TO_GAMMA);
    }

    public static function MoveTowards($current, $target , $maxDelta) {
        if ($maxDelta > 0) {
            if ($target < $current && $current - $maxDelta < $target) return $target;
            else if ($target > $current && $current + $maxDelta > $target) return $target;
        }
        if ($current > $target) {
            return $current - $maxDelta;
        }
        return $current + $maxDelta;
    }

    public static function PerlinNoise($x, $y) {
        return (new SimplexNoise(self::RANDOM_SEED))->Noise($x, $y);
    }

    public static function PingPong($t, $length) {
        if ($t < 0) $t = -$t;
        $mod = $t % $length;
        if (ceil($t / $length) % 2 === 0) {
            return ($mod === 0) ? 0 : $length - ($mod);
        }
        return ($mod === 0) ? $length : $mod;
    }

    public static function Repeat($t, $length) {
        if ($t > 0) return $t % $length;
        return $length + ($t % $length);
    }

    public static function Sign($f) {
        return ($f >= 0) ? 1 : -1;
    }


}
