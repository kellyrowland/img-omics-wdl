workflow crt {
  call run { }
}

task run {
  command {
	cmd="ls"
	cmd="$cmd * 1> /dev/null"
	$cmd
  }
}

