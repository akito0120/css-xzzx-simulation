from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn
from rich.console import Console
import time

console = Console()

with Progress(
    TextColumn("{task.description}"),
    BarColumn(bar_width=30),
    TextColumn("[green]{task.completed}/{task.total}"),
    console=console,
) as progress:

    projects = ["Project A", "Project B"]
    all_tasks = progress.add_task("[bold white] Simulation", total=len(projects))

    for project in projects:
        modules = ["Module 1", "Module 2", "Module 3"]
        proj_task = progress.add_task(f"[yellow]    {project}", total=len(modules))

        for module in modules:
            steps = 4
            mod_task = progress.add_task(f"[cyan]       {module}", total=steps)

            for _ in range(steps):
                time.sleep(0.2)
                progress.advance(mod_task)

            progress.advance(proj_task)
            progress.remove_task(mod_task)

        progress.advance(all_tasks)
        progress.remove_task(proj_task)
