use ExtUtils::MakeMaker;
WriteMakefile(
	      'NAME'      => 'JCMT::DAS',
	      'VERSION'   => '0.1',
	      'OBJECT'    => 'DAS.o dasmerge.o arrays.o',
	      'DEFINE'    => '',
    );
