use Random;
use Time;

proc monte_carlo_method(N: int, seed: int) : int {
	var randStreamSeeded = new owned RandomStream(real, seed);
	var count: int = 0; // number of points in the 1st quadrant of unit circle
	for i in 0..(N-1) {
		var x = randStreamSeeded.getNext();
		var y = randStreamSeeded.getNext();
		var z: real = x * x + y * y;
		if z <= 1 {
			count += 1;
		}
	}

	return count;
}

proc main(args: [] string){
	var num_threads: int = here.maxTaskPar;

	var global_count: int = 0; // number of points in the 1st quadrant of unit circle
	var N: int = 10000;
	if args.size >= 1 {
		N = args(1): int;
	}

	var timer = new Timer();
	timer.start();

	forall thread_id in 1..(num_threads) with (+ reduce global_count) do{
		global_count += monte_carlo_method(N/num_threads, thread_id);
	}

	var pi: real = global_count / N:real * 4;

	timer.stop();
	writeln("calculated value of pi:",pi);
	writeln("The process took ",timer.elapsed(TimeUnits.seconds)," seconds to finish");
	writeln("maxTaskPar:", num_threads);

	return;
}
