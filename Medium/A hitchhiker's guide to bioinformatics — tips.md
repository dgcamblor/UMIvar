---
subtitle: I've come across 1000 errors in my pipelines, here's what i learnt... or something like that would write your average clickbaiter
writingStatus: planning
---

Bioformatic tips I wish I learnt sooner?

# Planning

Title:

- Bioinformatics piques your interest? Here is advice for a biologist in need

First of all, start with a paragraph that expresses how bioinformatics is a great field, and how interesting it is specially with the recent shift to data generation, omics, etc. Make yourself clear that these advice will be mainly for biologists. 


- **Get yourself a comfortable environment.** Get a good IDE. Get yourself confortable in an OS (be it Windows or Ubuntu). If you use VIM, you may not be able to leave.

- **Learn the ways of scripting.** Scripts are the greatest power you have as a bioinformatician, and one you should embrace. They are not just for automating tasks, but also for documenting your work and making it reproducible.
- It may be straightforward to some, but the power of using computers and scripts allows us to automate pretty much any analysis we have in mind. Even if you are thinking that you are going to execute a pipeline just once, you should aim to automate things as if they were going to be repeated - chances are you will have to repeat them. Not only it reduces times, but it also reduces tedious and the mistakes it come with it.
- Coupling with a CONF file.

- **Adopt good programming practices beyond learning just to code.** Rembember to document your code, to write according to style guides and to use version control. This is both for your own sake and for the sake of others. We are not the same meme.

- **Keep a journal close to you.** It is a good idea to keep a journal of the commands you use, the problems you encounter and the solutions you find. This will help you to remember things and to learn from your mistakes. It is also a good idea to keep a cheatsheet of the commands you use the most. I personally use Obsidian for this purpose. But a plain notebook is just suficient.
- **Develop good debugging skills.** (?)

- Build up a good organizational folder system.

- **Sorround yourself with good people.** A goos balance between learning from others and learning by yourself is key to success. You should aim to be part of a community of bioinformaticians, either in person or online. This will help you to learn from others and to get help when you need it. It is also a good idea to have a mentor, someone who can guide you through the process of learning bioinformatics.

- **Get up to date with the latest technologies.** Last advice. Bioinformatics is an ever-evolving field (statistics, informatics, biology, etc.). Perhaps this is now more true than ever with the recent surge in AI-powered tools. Of course, no need to say that the purely "bio" part in bioinformatics is also constantly changing, and should not be neglected. Like the legendary Bruce Lee said, "Be water, my friend". Be flexible and prepere to adapt for the new ways of programming and doing bioinformatics that are coming (rather sooner than later).

I personally prefer inline code completers more than good ol' chatgpt, mainly because. But do take your time and see what works best for you. but critiquely

    - Good ol' ChatGPT
    - GitHub Copilot (free for students)

# Writing

Bioinformatics combines the best of a myriad of disciplines. I guess it's not too hard to figure out the first two —biology and informatics— but statistics, mathematics, etc. are also players in this game. The recent shift to data generation and omics has made bioinformatics a field more relevant than ever, and it is no wonder that many biologists are now interested in learning about it.

My career path was always clear to me - I wanted to become a biologist. But I grew practically attached to a computer, which made me also like that field. That is why when I learned about bioinformatics, I sort of found my Ikigai. Along my career path, I've made a lot of mistakes, but learned an equivalent amount of insights. If you are a biologist who is curious about, in need of or comitted to learning bioinformatics, I want to share you with you my knowledge!

## Get yourself a comfortable (computing) environment

You will probably know (if not, I am telling you now) that most bioinformatic tools are developed for Unix-like systems (say, Linux or MacOS), so you should get confortable with them. If you are not a MacOS user, the easiest way to get into Linux is by installing Ubuntu: it is just easy to do, and you get a greatly accessible environment to play with. However, if you are too in love with Windows (or you are a hardcore PC gamer), don't panic! There are plenty of options for you nowadays. You can use the WSL (Windows Subsystem for Linux) or set up a virtual machine (for which ther are plenty of tutorials online, I would personally recommend you to use VirtualBox).

Once tha matter of the OS is settled, you should get yourself a good code editor. Editors like Vim, Nano or Emacs can be great, but the steep learning curve is bound to make you give up, so I would only recommend you to use the first two for quick edits. A great beginner-friendly editor that is equally powerful and with a great community behind it is Visual Studio Code. It is available for all major OSs, and it has a vast ecosystem of plugins that make it highly customizable. And if you are working in a remote server, you can easily use it to connect to it and explore and edit your files. Although VSCode helps in writing in any language, for writing in R, RStudio is a more recommendable option, as it comes with a lot of features that are specific to the language (you can get them with VSCode, but needs a bit of tweaking).

Ah, let me not forget about mentioning the good ol' debate about the best language for bioinformatics. I will, however, be brief by mentioning that Bash, Python and R are specifically good for scripting, general purpose programming and statistics and data analysis, respectively. You do you, but it is probable that you will end up using all of them to some extent. If learning a language for the first time, Python is the easiest to 

## Learn the ways of scripting and follow good programming practices

Scripts are the greatest power you have as a bioinformatician, and one you should embrace. They not only offer you the ability to automate a set of steps, but they allow you to: 1) a record of each command you have launched, with the potential for debugging that so implies; and 2) offer a reproducible (or at least partially) way of reproducing a file processing or analysis you have conducted. Even if you think you will perform a set of steps just once, you should write them in a script. And chances are, you are going to perform them more times.

The documentation of the code you write helps both you and those who read it. In fact, documentation is one of 

## Have a journal close to you

## Keep up to date with the latest technologies