version 1.0
task t {
  input { String img }
  command <<< echo ok; uname -m >>>
  runtime { docker: img  maxRetries: 0  preemptible: 0 }
}
workflow w {
  input { String img }
  call t { input: img = img }
}
