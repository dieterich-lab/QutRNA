rule create_info_config:
  output:
    "info/config.yaml"
  run:
    with open(output[0], "w") as out:
      yaml.dump(config, out)


if workflow.configfiles:
  rule create_info_configfiles:
    input: workflow.configfiles
    output: [os.path.join("info/configfiles", f"{i}_{os.path.basename(cf)}") for i, cf in enumerate(workflow.configfiles, start=1)]
    run:
      for src, dst in zip(input, output):
        shutil.copyfile(src, dst)


if workflow.pepfile:
  rule create_info_pep:
    input: workflow.pepfile
    output:
      "info/pep.yaml"
    run:
      with open(output[0], "w") as out:
        yaml.dump(pep.to_dict(), out)

  rule create_info_pepfile:
    input: workflow.pepfile
    output:
      "info/pepfile.yaml"
    run:
      dst = os.path.basename(input[0])
      shutil.copyfile(src, dst)


rule create_info_version:
  output:
    "info/version.txt"
  run:
    qutrna_version = "QutRNA2: " + VERSION
    if shutil.which("git") and getattr(workflow, "_main_snakefile"):
      snakedir = os.path.dirname(workflow._main_snakefile)
      result = subprocess.run(["git", "status"], stdout=subprocess.PIPE, cwd=snakedir)
      qutrna_version += "\n" + result.stdout.decode("utf-8")
      result = subprocess.run(["git", "log", "-1"], stdout=subprocess.PIPE, cwd=snakedir)
      qutrna_version += "\n" + result.stdout.decode("utf-8")

    with open(output[0], "w") as out:
      out.write(qutrna_version)


rule create_info_opts:
  output:
    "info/opts.txt"
  run:
    with open(output[0], "w") as out:
      out.write(" ".join(sys.argv))
